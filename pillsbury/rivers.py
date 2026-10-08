"""Water depth in the rivers across the drained lakebed, for three measured flows.

Flows: USGS 11473900, Middle Fork Eel near Dos Rios (745 sq mi, no dam), scaled by drainage
area to every channel cell. The gauge below Scott Dam measures released water, not the
river; the Middle Fork is the nearest undammed record. Its median annual peak per square
mile (46 cfs) is close to the Scott Dam gauge's (41), which supports the scaling.

Depth: HAND with a synthetic rating curve per reach, the method behind NOAA's National
Water Model flood maps. HAND (height above nearest drainage) is each cell's height above
the channel cell it drains to. For a reach, filling every cell below stage h gives a
flow area and wetted perimeter, so Manning's equation gives Q(h). Inverting at the reach's
flow gives its stage, and depth = stage - HAND.
"""
import json
import numpy as np
import datetime as dt
from osgeo import gdal

gdal.UseExceptions()
MF_SQMI = 745
N_MANNING = 0.05   # bare sand and gravel with some roughness; NOAA uses 0.06 by default
MIN_SLOPE = 0.0005
CELL = 10.0
MIN_DEPTH_M = 0.02
CFS = 0.0283168


def rdb(path, date_col, value_col):
    out = {}
    for line in open(path):
        if line.startswith(("#", "agency", "5s")):
            continue
        p = line.rstrip("\n").split("\t")
        try:
            out[p[date_col]] = float(p[value_col])
        except (ValueError, IndexError):
            pass
    return out


def scenarios():
    daily = {dt.date.fromisoformat(k): v for k, v in rdb("data/dv_11473900.rdb", 2, 3).items()}
    peaks = rdb("data/pk_11473900.rdb", 2, 4)
    return {
        "summer": ("Late summer", np.median([q for d, q in daily.items() if d.month == 8])),
        "winter": ("Typical winter day", np.median([q for d, q in daily.items() if d.month in (1, 2)])),
        "flood": ("2-year flood peak", np.median(list(peaks.values()))),
    }


def read(path):
    ds = gdal.Open(path)
    a = ds.ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return a, ds


z, ref = read("data/dem10.tif")
acc, _ = read("data/acc10.tif")
drn, _ = read("data/dir10.tif")
seg, _ = read("data/streams.tif")
def on_grid(path, alg):
    g = ref.GetGeoTransform()
    return gdal.Warp("", path, format="MEM", outputBounds=(g[0], g[3] + g[5] * ref.RasterYSize, g[0] + g[1] * ref.RasterXSize, g[3]),
                     width=ref.RasterXSize, height=ref.RasterYSize, resampleAlg=alg).ReadAsArray() > 0


basin = on_grid("data/basin.tif", "mode")
notch = on_grid("data/notch.tif", "max")
H, W = z.shape
n = H * W

# downstream neighbour of every cell (r.watershed codes 1..8 = NE, N, NW, W, SW, S, SE, E)
dr = np.array([0, -1, -1, -1, 0, 1, 1, 1, 0])
dc = np.array([0, 1, 0, -1, -1, -1, 0, 1, 1])
code = np.nan_to_num(drn, nan=0).astype(int)
rr, cc = np.indices(z.shape)
r2, c2 = rr + dr[np.abs(code)], cc + dc[np.abs(code)]
ok = (code > 0) & (r2 >= 0) & (r2 < H) & (c2 >= 0) & (c2 < W)
down = np.where(ok, r2 * W + c2, -1).ravel()

# nearest channel cell downstream, by pointer doubling: channel cells point at themselves,
# so each jump doubles the distance followed until every path has stopped on one
is_stream = ~np.isnan(seg.ravel())
target = np.where(is_stream, np.arange(n), down)
while True:
    nxt = np.where(target >= 0, target[np.maximum(target, 0)], -1)
    if np.array_equal(nxt, target):
        break
    target = nxt
reached = (target >= 0) & is_stream[np.where(target >= 0, target, 0)]
zf = z.ravel()
hand = np.where(reached, zf - zf[np.where(reached, target, 0)], np.nan)
reach = np.where(reached, seg.ravel()[np.where(reached, target, 0)], np.nan)

gy, gx = np.gradient(z, CELL)
slope_f = np.sqrt(1 + gx ** 2 + gy ** 2).ravel()
area_m2 = np.abs(acc.ravel()) * CELL ** 2

# reach length, slope and drainage area at its downstream end
ids = np.unique(seg[~np.isnan(seg)]).astype(int)
segf = seg.ravel()
reaches = {}
for i in ids:
    cells = np.nonzero(segf == i)[0]
    L = len(cells) * CELL * 1.13  # mean step of a D8 path
    zc = zf[cells]
    reaches[i] = dict(L=L, S=max((np.nanmax(zc) - np.nanmin(zc)) / L, MIN_SLOPE),
                      A=area_m2[cells].max())
# fill-stage tables: per reach, its cells' HAND sorted, with running sums
order = np.argsort(np.nan_to_num(reach, nan=-1) * 1e6 + np.nan_to_num(hand, nan=1e5), kind="stable")
groups = {}
r_sorted = reach[order]
valid = ~np.isnan(r_sorted)
starts = np.r_[0, np.nonzero(np.diff(r_sorted[valid]))[0] + 1]
ov = order[valid]
for s, e in zip(starts, np.r_[starts[1:], valid.sum()]):
    cells = ov[s:e]
    h = hand[cells]
    keep = h < 25
    cells, h = cells[keep], h[keep]
    groups[int(reach[cells[0]])] = (cells, h, np.cumsum(h), np.cumsum(slope_f[cells]))


def stage(i, q):
    r = reaches[i]
    cells, h, csum, cf = groups[i]
    stages = np.arange(0.01, 25, 0.01)
    k = np.searchsorted(h, stages)
    ok = k > 0
    V = np.where(ok, k * stages - np.where(ok, csum[np.maximum(k - 1, 0)], 0), 0) * CELL ** 2
    P = np.where(ok, cf[np.maximum(k - 1, 0)], 1) * CELL ** 2 / r["L"]
    A = V / r["L"]
    Q = A * (A / P) ** (2 / 3) * np.sqrt(r["S"]) / N_MANNING
    j = np.searchsorted(Q, q)
    return stages[min(j, len(stages) - 1)]


grass = dict(l.strip().split("=") for l in open("data/rivers_grass.txt"))
upper_eel_sqmi = float(grass["upper_eel_km2"]) / 2.58999
show = (basin | notch).ravel()
dam_acc = area_m2[notch.ravel()].max() / 2.58999e6  # drainage area where the river leaves the notch

lake = gdal.Open("web/lakebed.f32")
lgt, LW, LH = lake.GetGeoTransform(), lake.RasterXSize, lake.RasterYSize
bounds = (lgt[0], lgt[3] + lgt[5] * LH, lgt[0] + lgt[1] * LW, lgt[3])


def to_web(arr, path, resample):
    mem = gdal.GetDriverByName("MEM").Create("", W, H, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(ref.GetGeoTransform())
    mem.SetProjection(ref.GetProjection())
    mem.GetRasterBand(1).SetNoDataValue(-9999)
    mem.GetRasterBand(1).WriteArray(np.nan_to_num(arr, nan=-9999).astype(np.float32))
    return gdal.Warp(path, mem, dstSRS="EPSG:4326", outputBounds=bounds, width=LW, height=LH,
                     resampleAlg=resample, dstNodata=-9999, format="MEM" if path == "" else "ENVI")


meta = json.load(open("web/meta.json"))
meta["rivers"] = {}
for key, (label, cfs_mf) in scenarios().items():
    q_per_m2 = cfs_mf * CFS / (MF_SQMI * 2.58999e6)
    depth = np.full(n, np.nan)
    for i in groups:
        hs = stage(i, q_per_m2 * reaches[i]["A"])
        cells, h, _, _ = groups[i]
        wet = h < hs
        depth[cells[wet]] = hs - h[wet]
    depth[~show | (depth < MIN_DEPTH_M)] = np.nan
    d = depth.reshape(H, W)
    web = to_web(d, "", "bilinear").ReadAsArray()
    web[web < MIN_DEPTH_M] = np.nan
    ft = web / 0.3048006
    # pale at a few inches, mid blue at 5 ft, navy from 20 ft
    t1, t2 = np.clip(ft / 5, 0, 1), np.clip((ft - 5) / 15, 0, 1)
    rgb = [np.where(ft < 5, a + (b - a) * t1, b + (c - b) * t2) for a, b, c in ((170, 40, 10), (220, 120, 40), (250, 200, 110))]
    rgba = np.stack(rgb + [np.where(np.isnan(ft), 0, 255)]).astype(np.uint8)
    for k in range(3):
        rgba[k][np.isnan(ft)] = 0
    png = gdal.GetDriverByName("MEM").Create("", LW, LH, 4, gdal.GDT_Byte)
    for k in range(4):
        png.GetRasterBand(k + 1).WriteArray(rgba[k])
    gdal.GetDriverByName("PNG").CreateCopy(f"web/river_{key}.png", png)
    wet_acres = np.isfinite(d).sum() * CELL ** 2 / 4046.86
    cfs_dam = cfs_mf * upper_eel_sqmi / MF_SQMI
    meta["rivers"][key] = dict(label=label, cfs_at_dam=round(float(cfs_dam)), wet_acres=round(float(wet_acres)),
                               max_depth_ft=round(float(np.nanmax(d) / 0.3048006), 1))
    print(f"{label}: {cfs_dam:.0f} cfs at the dam; {wet_acres:.0f} acres wet in the basin; deepest {np.nanmax(d) / 0.3048006:.1f} ft")

meta["upper_eel_sqmi"] = round(upper_eel_sqmi)
meta["drainage_at_dam_sqmi_5m"] = round(dam_acc)
json.dump(meta, open("web/meta.json", "w"), indent=1)
print(f"upper Eel above the dam: {upper_eel_sqmi:.0f} sq mi (USGS: 290); routed through the notch at 10 m: {dam_acc:.0f} sq mi")

# the terrain the viewer draws: the drained DEM, notch included
gdal.Warp("web/lakebed.f32", "data/drained_navd88_m.tif", dstSRS="EPSG:4326", outputBounds=bounds, width=LW, height=LH,
          resampleAlg="bilinear", dstNodata=-9999, format="ENVI", outputType=gdal.GDT_Float32)
