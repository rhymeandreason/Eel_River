"""PG&E 2023 Lake Pillsbury survey PDF -> elevation-labelled contour lines.

The survey's Sheet 1 is an AutoCAD export: 5 ft contours as vector paths on layer
"Surface 2023" (every 5th, the 25 ft majors, in a darker colour) over a colour-ramped
bathymetry raster on layer "2023 Bathy Image". There are no elevation labels and no
geospatial tags, so this script recovers both:

  position   scale from the sheet (1 pt = 27.78 ft), then shift/rotate fitted so the
             contours above the lidar's water line run level on the 3DEP lidar
  datum      PG&E datum = NAVD88 + OFF ft, fitted from those same dry contours, which
             must land on multiples of 5 (and the majors on multiples of 25)
  elevation  dry contours: read off the lidar. Submerged: the ramp colour under the
             line, through HUE_LADDER below

Writes data/contours.gpkg, data/bathy_hue.tif (ramp read as elevation, for QA),
data/footprint.tif (where the survey has data) and data/extract.json.
"""
import json, sys
from collections import defaultdict
import numpy as np
import pymupdf
from osgeo import gdal, ogr, osr
from scipy import ndimage, optimize, signal

gdal.UseExceptions()
PDF, LIDAR = "data/bathy2023.pdf", "data/lidar_ft.tif"
PAGE = 12          # Sheet 1, 0-based
FT_PER_PT = 27.7778  # viewport /Measure and the 0-4000 ft scale bar agree
DPI = 300
US_FT = 1200 / 3937

# Ramp hue (deg) at each contour level, PG&E ft. Read from the clusters of contour hues:
# 25 ft majors fall at 66, 190 and 221 with exactly four minor clusters between each pair,
# and the 1890 cluster (hue 20) checks against the lidar where both exist. The blue end
# is compressed (about 4 deg per 5 ft below 1850), so labels there are the least certain.
HUE_LADDER = [(13, 1895), (20, 1890), (33, 1885), (50, 1880), (66, 1875), (80, 1870),
              (97, 1865), (132, 1860), (165, 1855), (190, 1850), (198, 1845), (204, 1840),
              (210, 1835), (215, 1830), (221, 1825), (234, 1820), (240, 1815)]


def page_contours(page):
    """Contour fragments from the Surface layer, chained end to end, in display points
    (north up, y down)."""
    rot = page.rotation_matrix
    frags = []
    for d in page.get_drawings():
        if d.get("layer") != "Surface 2023":
            continue
        major = d["color"][0] < 0.4
        pts = []
        for it in d["items"]:
            if it[0] not in ("l", "c"):
                continue
            a, b = (it[1], it[2]) if it[0] == "l" else (it[1], it[4])
            a, b = a * rot, b * rot
            if pts and abs(pts[-1][0] - a.x) < 1e-3 and abs(pts[-1][1] - a.y) < 1e-3:
                pts.append((b.x, b.y))
            else:
                if len(pts) > 1:
                    frags.append((major, np.array(pts)))
                pts = [(a.x, a.y), (b.x, b.y)]
        if len(pts) > 1:
            frags.append((major, np.array(pts)))

    key = lambda q: (round(q[0], 2), round(q[1], 2))
    ends = defaultdict(list)
    for i, (m, p) in enumerate(frags):
        ends[(m, key(p[0]))].append(i)
        ends[(m, key(p[-1]))].append(i)
    used = np.zeros(len(frags), bool)
    chains = []
    for i in range(len(frags)):
        if used[i]:
            continue
        used[i] = True
        m, p = frags[i]
        line = [tuple(q) for q in p]
        for _ in (0, 1):
            while True:
                k = (m, key(line[-1]))
                nxt = [j for j in ends[k] if not used[j]]
                if not nxt:
                    break
                used[nxt[0]] = True
                q = frags[nxt[0]][1]
                q = q if key(q[0]) == k[1] else q[::-1]
                line += [tuple(v) for v in q[1:]]
            line.reverse()
        chains.append((m, np.array(line)))
    return chains


def ramp_image(doc, page):
    """Render only the bathymetry raster layer; return hue (deg) and saturation."""
    for c in doc.layer_ui_configs():
        if c["type"] == "checkbox":
            doc.set_layer_ui_config(c["number"], 0 if c["text"] == "2023 Bathy Image" else 2)
    pm = page.get_pixmap(dpi=DPI)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[..., :3] / 255.0
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    mx, mn = a.max(-1), a.min(-1)
    d = np.maximum(mx - mn, 1e-6)
    hue = np.where(mx == r, ((g - b) / d) % 6, np.where(mx == g, (b - r) / d + 2, (r - g) / d + 4)) * 60
    sat = (mx - mn) / np.maximum(mx, 1e-6)
    return hue.astype(np.float32), sat.astype(np.float32)


class Lidar:
    def __init__(self, path):
        ds = gdal.Open(path)
        self.z = ds.ReadAsArray().astype(float)
        self.z[self.z < -1000] = np.nan
        self.gt = ds.GetGeoTransform()
        v, n = np.unique(np.round(self.z[~np.isnan(self.z)], 1), return_counts=True)
        self.water_level = float(v[n.argmax()])  # the flat lake surface is the commonest value
        gy, gx = np.gradient(self.z, self.gt[1])
        self.slope = np.hypot(gx, gy)

    def rc(self, X, Y):
        return [(self.gt[3] - Y) / -self.gt[5] - 0.5, (X - self.gt[0]) / self.gt[1] - 0.5]

    def sample(self, X, Y, grid=None):
        return ndimage.map_coordinates(self.z if grid is None else grid, self.rc(X, Y), order=1, cval=np.nan)


def to_world(c, g):
    th = np.deg2rad(g["rot"])
    s = FT_PER_PT * g["k"]
    x, y = c[:, 0] * s, -c[:, 1] * s
    return g["X0"] + x * np.cos(th) - y * np.sin(th), g["Y0"] + x * np.sin(th) + y * np.cos(th)


def georeference(chains, lid):
    """Coarse: slide the contours' footprint over the lidar's lake. Fine: make the
    contours above the water line as level as possible on the lidar."""
    cell = lid.gt[1]
    water = ndimage.binary_opening(np.abs(lid.z - lid.water_level) < 0.15, iterations=2)
    lab, _ = ndimage.label(water)
    water = lab == np.bincount(lab.ravel())[1:].argmax() + 1

    pts = np.vstack([c for _, c in chains])
    F = np.zeros((int(800 * FT_PER_PT / cell) + 1, int(1224 * FT_PER_PT / cell) + 1), bool)
    F[(pts[:, 1] * FT_PER_PT / cell).astype(int), (pts[:, 0] * FT_PER_PT / cell).astype(int)] = True
    F = ndimage.binary_fill_holes(ndimage.binary_closing(F, iterations=15))
    rows, cols = np.nonzero(F)
    T = F[rows.min():rows.max() + 1, cols.min():cols.max() + 1].astype(float)
    corr = signal.fftconvolve(np.where(water, 1.0, -1.0), T[::-1, ::-1], mode="valid")
    iy, ix = np.unravel_index(corr.argmax(), corr.shape)
    g0 = dict(X0=lid.gt[0] + (ix - cols.min()) * cell, Y0=lid.gt[3] - (iy - rows.min()) * cell, rot=0.0, k=1.0)

    sel = []
    for _, c in chains:
        d = np.r_[0, np.cumsum(np.hypot(*np.diff(c, axis=0).T))]
        if d[-1] >= 15:
            t = np.arange(0, d[-1], 0.5)
            sel.append(np.c_[np.interp(t, d, c[:, 0]), np.interp(t, d, c[:, 1])])

    def cost(p):
        g = dict(X0=g0["X0"] + p[0], Y0=g0["Y0"] + p[1], rot=p[2], k=1 + p[3])
        sds = []
        for c in sel:
            z = lid.sample(*to_world(c, g))
            dry = z > lid.water_level + 1
            if dry.mean() > 0.9 and dry.sum() > 30:
                sds.append(np.std(z[dry]))
        return np.median(sds) if len(sds) > 20 else 99.0

    best = min(((cost([dx, dy, 0, 0]), [dx, dy, 0, 0]) for dx in range(-400, 401, 50) for dy in range(-400, 401, 50)))
    simplex = np.array(best[1], float) + np.array([[0, 0, 0, 0], [30, 0, 0, 0], [0, 30, 0, 0], [0, 0, .3, 0], [0, 0, 0, .003]])
    r = optimize.minimize(cost, best[1], method="Nelder-Mead",
                          options=dict(initial_simplex=simplex, xatol=0.5, fatol=1e-3, maxiter=400))
    p = r.x
    return dict(X0=g0["X0"] + p[0], Y0=g0["Y0"] + p[1], rot=float(p[2]), k=float(1 + p[3]), level_sd_ft=float(r.fun))


def hue_to_pge(h):
    hs, es = zip(*HUE_LADDER)
    return np.interp(h, hs, es, left=np.nan, right=np.nan)


def snap(e, major):
    if major:
        return 25 * round(e / 25)
    v = 5 * round(e / 5)
    if v % 25 == 0:  # a minor line is never on a major level
        v += 5 if e > v else -5
    return v


def main():
    doc = pymupdf.open(PDF)
    page = doc[PAGE]
    chains = page_contours(page)
    lid = Lidar(LIDAR)
    print(f"{len(chains)} contour lines; lidar water surface {lid.water_level} ft NAVD88", file=sys.stderr)

    g = georeference(chains, lid)
    print("georeference", g, file=sys.stderr)

    # per-line lidar reading where the line is on dry, gentle ground
    lines = []
    for m, c in chains:
        X, Y = to_world(c, g)
        z, sl = lid.sample(X, Y), lid.sample(X, Y, lid.slope)
        dry = (z > lid.water_level + 1) & (sl < 0.25)
        zl = np.median(z[dry]) if dry.sum() > 15 and dry.mean() > 0.9 and np.std(z[dry]) < 1.5 else np.nan
        lines.append(dict(major=bool(m), page=c, X=X, Y=Y, z_lidar=zl))

    # datum offset: dry lines must sit on 5 ft multiples
    zl = np.array([l["z_lidar"] for l in lines if not np.isnan(l["z_lidar"])])
    offs = np.arange(76, 82, 0.05)
    miss = [np.mean(np.abs((zl + o + 2.5) % 5 - 2.5)) for o in offs]
    OFF = float(offs[int(np.argmin(miss))])
    maj = np.array([l["z_lidar"] + OFF for l in lines if l["major"] and not np.isnan(l["z_lidar"])])
    print(f"PG&E = NAVD88 + {OFF:.2f} ft; dry lines {len(zl)}, mean miss {min(miss):.2f} ft; "
          f"dry majors read {np.round(np.sort(maj)).astype(int).tolist()}", file=sys.stderr)

    hue, sat = ramp_image(doc, page)
    k = DPI / 72
    labelled, agree, compared = [], 0, 0
    for l in lines:
        px = np.clip((l["page"] * k).astype(int), 0, [hue.shape[1] - 1, hue.shape[0] - 1])
        hh, ss = hue[px[:, 1], px[:, 0]], sat[px[:, 1], px[:, 0]]
        e_hue = hue_to_pge(np.median(hh[ss > 0.45])) if (ss > 0.45).mean() > 0.6 else np.nan
        e_lid = l["z_lidar"] + OFF if not np.isnan(l["z_lidar"]) else np.nan
        if not np.isnan(e_hue) and not np.isnan(e_lid):
            compared += 1
            agree += snap(e_hue, l["major"]) == snap(e_lid, l["major"])
        if not np.isnan(e_lid):
            l.update(pge=snap(e_lid, l["major"]), src="lidar")
        elif not np.isnan(e_hue):
            l.update(pge=snap(e_hue, l["major"]), src="ramp")
        else:
            continue
        labelled.append(l)
    print(f"labelled {len(labelled)} of {len(lines)}; ramp vs lidar agree on {agree}/{compared}", file=sys.stderr)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(2226)
    drv = ogr.GetDriverByName("GPKG")
    out = drv.CreateDataSource("data/contours.gpkg")
    lyr = out.CreateLayer("contours", srs, ogr.wkbLineString)
    for name, t in (("pge_ft", ogr.OFTReal), ("navd88_m", ogr.OFTReal), ("major", ogr.OFTInteger), ("src", ogr.OFTString)):
        lyr.CreateField(ogr.FieldDefn(name, t))
    for l in labelled:
        f = ogr.Feature(lyr.GetLayerDefn())
        f["pge_ft"], f["navd88_m"] = l["pge"], (l["pge"] - OFF) * US_FT
        f["major"], f["src"] = int(l["major"]), l["src"]
        geom = ogr.Geometry(ogr.wkbLineString)
        for x, y in zip(l["X"], l["Y"]):
            geom.AddPoint_2D(float(x), float(y))
        f.SetGeometry(geom)
        lyr.CreateFeature(f)
    out = None

    # the ramp itself, georeferenced, as NAVD88 m (QA) and as the survey footprint
    th, s = np.deg2rad(g["rot"]), FT_PER_PT * g["k"] / k
    gt = (g["X0"], s * np.cos(th), s * np.sin(th), g["Y0"], s * np.sin(th), -s * np.cos(th))
    e = (hue_to_pge(hue) - OFF) * US_FT
    e[sat < 0.45] = np.nan
    for path, arr, typ, nd in (("data/bathy_hue.tif", e, gdal.GDT_Float32, -9999),
                               ("data/footprint.tif", (sat >= 0.45).astype(np.uint8), gdal.GDT_Byte, 0)):
        ds = gdal.GetDriverByName("GTiff").Create(path, arr.shape[1], arr.shape[0], 1, typ, ["COMPRESS=DEFLATE"])
        ds.SetGeoTransform(gt)
        ds.SetProjection(srs.ExportToWkt())
        b = ds.GetRasterBand(1)
        b.SetNoDataValue(nd)
        b.WriteArray(np.nan_to_num(arr, nan=nd))
        ds = None

    json.dump(dict(georef=g, pge_minus_navd88_ft=OFF, lidar_water_ft_navd88=lid.water_level,
                   lines=len(lines), labelled=len(labelled), ramp_lidar_agree=[int(agree), compared],
                   by_source={s_: sum(l["src"] == s_ for l in labelled) for s_ in ("lidar", "ramp")}),
              open("data/extract.json", "w"), indent=1)


if __name__ == "__main__":
    main()
