"""Take Scott Dam out of the lakebed DEM.

The bare-earth lidar still holds a ridge about 1890 ft high at the dam site (what is left
of the dam and its footings after ground classification), so a drained basin would keep a
pond 75 ft deep. Removal takes the dam down to the river, so this cuts a notch from the
lowest lake floor near the dam to the lowest river bed below it, on a straight grade, and
only ever lowers the ground. The real removal design (and how much sediment is dug out at
the dam) may differ.
"""
import json
import numpy as np
from osgeo import gdal

gdal.UseExceptions()
NOTCH_HALF_WIDTH_M = 10
SEARCH_M = 500

dam = json.load(open("data/dam.json"))
ds = gdal.Open("data/lakebed_navd88_m.tif")
z = ds.ReadAsArray().astype(float)
nd = ds.GetRasterBand(1).GetNoDataValue()
gt = ds.GetGeoTransform()
cell = gt[1]
basin = gdal.Open("data/basin.tif").ReadAsArray() == 1

rows, cols = np.indices(z.shape)
X, Y = gt[0] + (cols + 0.5) * cell, gt[3] - (rows + 0.5) * cell
near = np.hypot(X - dam["utm_x"], Y - dam["utm_y"]) < SEARCH_M
a = np.unravel_index(np.argmin(np.where(near & basin, z, np.inf)), z.shape)
b = np.unravel_index(np.argmin(np.where(near & ~basin & (z != nd), z, np.inf)), z.shape)

# distance along and across the line a -> b, for every cell
ax, ay, bx, by = X[a], Y[a], X[b], Y[b]
L = np.hypot(bx - ax, by - ay)
t = ((X - ax) * (bx - ax) + (Y - ay) * (by - ay)) / L ** 2
d = np.abs((X - ax) * (by - ay) - (Y - ay) * (bx - ax)) / L
notch = (t >= 0) & (t <= 1) & (d <= NOTCH_HALF_WIDTH_M)
grade = z[a] + (z[b] - z[a]) * np.clip(t, 0, 1)
out = np.where(notch, np.minimum(z, grade), z)

ft = lambda m: m / 0.3048006 + json.load(open("data/extract.json"))["pge_minus_navd88_ft"]
print(f"notch {L:.0f} m from lake floor {ft(z[a]):.0f} ft to river bed {ft(z[b]):.0f} ft PG&E; "
      f"highest ground cut: {ft(z[notch].max()):.0f} ft")
o = gdal.GetDriverByName("GTiff").CreateCopy("data/drained_navd88_m.tif", ds, options=["COMPRESS=DEFLATE"])
o.GetRasterBand(1).WriteArray(out.astype(np.float32))
o = None
m = gdal.GetDriverByName("GTiff").CreateCopy("data/notch.tif", gdal.Open("data/basin.tif"), options=["COMPRESS=DEFLATE"])
m.GetRasterBand(1).WriteArray(notch.astype(np.uint8))
m = None
