"""Gravelly Valley before the dam, from the 1916 Army survey, as two layers for the viewer.

Source: USGS 1:62,500 Hullville quadrangle, published 1922: a U.S. Army Corps of Engineers
tactical map, "Surveyed in 1916", before Scott Dam flooded the valley. USGS serves it as a
GeoPDF (polyconic on NAD27), which GDAL reads with the libgdal-pdf plugin.

1. Warp the sheet onto the viewer's lake grid at twice its resolution.
2. Shift it to fit. A reconnaissance survey is good to about 100 m. The shift that best
   lays the sheet's blue ink on streams derived from the lidar, outside the lake where the
   ground hasn't changed, is found by brute-force correlation and applied.
3. web/map1916.webp: the sheet, kept only inside the lake basin and where the party mapped.
4. web/river1916.webp: the blue ink alone (the rivers and creeks of 1916) in the basin.

Run after build.sh and rivers.sh (it needs web/basin.u8 and data/streams.tif):
    conda run -n grass python hist.py
"""
import json, os, urllib.request
import numpy as np
from osgeo import gdal
from PIL import Image
from scipy import ndimage, signal

gdal.UseExceptions()
URL = ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Maps/HistoricalTopo/PDF/CA/62500/"
       "CA_Hullville_465755_1922_62500_geo.pdf")
PDF = "data/hist/hullville_1922.pdf"
UP = 2          # output pixels per lake-grid cell
PAD = 80        # search margin, in output pixels (about 250 m)
INK = (24, 66, 128)    # the traced channels: dark river blue
HALO = (246, 242, 233)

if not os.path.exists(PDF):
    os.makedirs(os.path.dirname(PDF), exist_ok=True)
    urllib.request.urlretrieve(URL, PDF)

meta = json.load(open("web/meta.json"))
G = meta["grid"]
W, H = G["width"] * UP, G["height"] * UP
px, py = (G["east"] - G["west"]) / W, (G["north"] - G["south"]) / H
bounds = (G["west"] - PAD * px, G["south"] - PAD * py, G["east"] + PAD * px, G["north"] + PAD * py)
size = dict(width=W + 2 * PAD, height=H + 2 * PAD, outputBounds=bounds, dstSRS="EPSG:4326")

gdal.SetConfigOption("GDAL_PDF_DPI", "400")   # about 4 m per pixel on the sheet
sheet = gdal.Warp("", PDF, format="MEM", resampleAlg="cubic", **size).ReadAsArray()
streams = gdal.Warp("", "data/streams.tif", format="MEM", resampleAlg="max", **size).ReadAsArray()
r, g, b = sheet.astype(int)
blue = (b - r > 40) & (b > g)          # river ink, strict: for the fit
water = (b - r > 25) & (b >= g)        # also the pale fill between a wide river's banks
streams = (streams > 0) & (streams < 1e9)

basin = np.fromfile("web/basin.u8", np.uint8).reshape(G["height"], G["width"]) > 0
basin = np.pad(basin.repeat(UP, 0).repeat(UP, 1), PAD)
# where the 1916 party mapped: dense contour and vegetation ink, holes filled
surveyed = ndimage.uniform_filter(((r < 200) | (g < 200)).astype(float), 101) > 0.08
surveyed = ndimage.binary_fill_holes(ndimage.binary_closing(surveyed, iterations=20))

# 2. compare only outside the lake, where the 1916 streams and today's still share ground
cmp = surveyed & ~ndimage.binary_dilation(basin, iterations=20)
A = ndimage.gaussian_filter((blue & cmp).astype(float), 6)
B = ndimage.gaussian_filter((streams & cmp).astype(float), 6)
# c[k] = sum B[n + s] A[n] at lag s = k - (N - 1): the score of rolling A by s
c = signal.correlate(B, A, mode="full", method="fft")
c0 = np.array(A.shape) - 1
win = c[c0[0] - PAD + 1:c0[0] + PAD, c0[1] - PAD + 1:c0[1] + PAD]
dy, dx = np.array(np.unravel_index(win.argmax(), win.shape)) - PAD + 1
scores = {(dx, dy): win.max(), (0, 0): c[tuple(c0)]}
m_per_px = 111320 * np.array([px * np.cos(np.radians((G["north"] + G["south"]) / 2)), py])
shift_e, shift_n = float(dx * m_per_px[0]), float(-dy * m_per_px[1])
print(f"shift {shift_e:+.0f} m east, {shift_n:+.0f} m north; "
      f"fit {scores[dx, dy] / scores[0, 0]:.2f}x the unshifted sheet")

# 2b. The survey's error grows away from its control (in the east arm the 1916 river runs
# 100-150 m south of the valley), so fit again tile by tile, now inside the lake too where
# its narrow arms still hold the old channel, within LOCAL px of the global shift. Tiles
# with little ink or no clear peak say nothing; the rest blend into a smooth warp that
# relaxes back to the global shift away from them.
A = ndimage.gaussian_filter((blue & surveyed).astype(float), 6)
B = ndimage.gaussian_filter((streams & surveyed).astype(float), 6)
T, LOCAL = 360, 50
num, den = np.zeros((2, *A.shape)), np.zeros(A.shape)
for y in range(LOCAL + abs(dy), A.shape[0] - T - LOCAL - abs(dy), T // 2):
    for x in range(LOCAL + abs(dx), A.shape[1] - T - LOCAL - abs(dx), T // 2):
        b_ = B[y:y + T, x:x + T]
        a_ = A[y - dy - LOCAL:y - dy + T + LOCAL, x - dx - LOCAL:x - dx + T + LOCAL]
        if b_.sum() < 200 or a_.sum() < 200:
            continue
        out = signal.correlate(a_, b_, mode="valid", method="fft")   # out[i, j]: shift L - i, L - j
        if out.max() < 1.5 * np.median(out):
            continue
        i, j = np.unravel_index(out.argmax(), out.shape)
        cy, cx = y + T // 2, x + T // 2
        num[:, cy, cx] = LOCAL - i, LOCAL - j
        den[cy, cx] = 1
sig = T / 2
w = ndimage.gaussian_filter(den, sig)
field = np.array([ndimage.gaussian_filter(n, sig) for n in num]) / (w + 0.2 * w.max())
local_m = float(np.hypot(*(field[::-1] * m_per_px[:, None, None])).max())
print(f"local correction up to {local_m:.0f} m from {int(den.sum())} tiles")

# sample the padded sheet through the warp: output pixel n shows source n - shift
Y, X = np.mgrid[PAD:PAD + H, PAD:PAD + W].astype(np.float32)
src = [Y - dy - field[0, PAD:PAD + H, PAD:PAD + W], X - dx - field[1, PAD:PAD + H, PAD:PAD + W]]
warp = lambda a, order: ndimage.map_coordinates(a, src, order=order, mode="nearest")
sheet = np.dstack([warp(band, 1) for band in sheet]).astype(np.uint8)
water, surveyed = warp(water, 0), warp(surveyed, 0)
basin = basin[PAD:PAD + H, PAD:PAD + W]

# 3. the sheet, feathered out a little past the full-pool shore
reach = ndimage.binary_dilation(basin, iterations=12) & surveyed
alpha = ndimage.gaussian_filter(reach.astype(float), 4)
Image.fromarray(np.dstack([sheet, (alpha * 255).astype(np.uint8)])).save("web/map1916.webp", quality=82)

# 4. the blue ink: bridge the breaks where labels and grid lines cross a river, fill
# between a wide river's banks (but not round an island), drop specks, then draw it in one
# colour on a pale halo
ink = ndimage.binary_closing(water, iterations=3) & ndimage.binary_dilation(basin, iterations=20)
lab, n = ndimage.label(ndimage.binary_fill_holes(ink) & ~ink)
small = 1 + np.flatnonzero(ndimage.sum(np.ones_like(lab), lab, range(1, n + 1)) < 400)
ink |= np.isin(lab, small)
lab, n = ndimage.label(ink, np.ones((3, 3)))
ink = np.isin(lab, 1 + np.flatnonzero(ndimage.sum(ink, lab, range(1, n + 1)) >= 150))
halo = ndimage.binary_dilation(ink, iterations=2) & ~ink
out = np.zeros((H, W, 4), np.uint8)
out[halo] = (*HALO, 150)
out[ink] = (*INK, 255)
Image.fromarray(out).save("web/river1916.webp", lossless=True)

meta["survey1916"] = dict(
    source="USGS 1:62,500 Hullville quadrangle (1922), U.S. Army Corps of Engineers, surveyed 1916",
    shift_east_m=round(shift_e), shift_north_m=round(shift_n), local_correction_max_m=round(local_m),
)
json.dump(meta, open("web/meta.json", "w"), indent=1)
