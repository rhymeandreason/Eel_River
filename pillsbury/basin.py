"""Where the reservoir can hold water: the cells below full pool that connect to the
lake bottom, minus the canyon below Scott Dam.

The 3DEP lidar is bare earth, so Scott Dam is not in it, and water above about 1891 ft
would flow straight down the Eel. The canyon is told apart from the lakebed by depth:
it holds ground lower than anything the 2023 survey found in the lake.
"""
import json
import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()
FT = 1200 / 3937
DAM_PLUG_M = 150  # radius blocked around the dam site
off = json.load(open("data/extract.json"))["pge_minus_navd88_ft"]
ds = gdal.Open("data/lakebed_navd88_m.tif")
z = ds.ReadAsArray().astype(float)
z[z == ds.GetRasterBand(1).GetNoDataValue()] = np.nan
surveyed = gdal.Open("data/footprint_utm.tif").ReadAsArray() == 1
lake_floor = np.nanmin(np.where(surveyed, z, np.nan))
seed = np.unravel_index(np.nanargmin(np.where(surveyed, z, np.nan)), z.shape)


def flood(pge, blocked=None):
    below = z < (pge - off) * FT
    if blocked is not None:
        below &= ~blocked
    lab, _ = ndimage.label(below)
    return lab == lab[seed]


# The lake at the lidar's water line is still held by the canyon walls. Below it, the
# canyon is the ground under that level that the lake does not reach and that dips below
# the lake floor; the dam stood where the canyon comes closest to the lake.
contained = flood(1888)
low = z < (1888 - off) * FT
lab, _ = ndimage.label(low & ~contained)
canyon = np.isin(lab, np.unique(lab[(lab > 0) & (z < lake_floor)]))
dist = ndimage.distance_transform_edt(~contained)
dam = np.unravel_index(np.argmin(np.where(canyon, dist, np.inf)), z.shape)
cell = ds.GetGeoTransform()[1]
yy, xx = np.ogrid[:z.shape[0], :z.shape[1]]
plug = (np.hypot(yy - dam[0], xx - dam[1]) * cell < DAM_PLUG_M) & ~contained
basin = flood(1910, blocked=canyon | plug)
leak = np.nanmin(np.where(basin, z, np.nan)) < lake_floor - 0.01
x, y = ds.GetGeoTransform()[0] + (dam[1] + 0.5) * cell, ds.GetGeoTransform()[3] - (dam[0] + 0.5) * cell
print(f"lake floor {lake_floor / FT + off:.1f} ft PG&E; dam plug at UTM {x:.0f},{y:.0f}; "
      f"basin {basin.sum() * cell ** 2 / 4046.86:.0f} acres at 1910; leaks: {leak}")
assert not leak
json.dump(dict(utm_x=float(x), utm_y=float(y), plug_m=DAM_PLUG_M), open("data/dam.json", "w"))

out = gdal.GetDriverByName("GTiff").Create("data/basin.tif", z.shape[1], z.shape[0], 1, gdal.GDT_Byte, ["COMPRESS=DEFLATE"])
out.SetGeoTransform(ds.GetGeoTransform())
out.SetProjection(ds.GetProjection())
out.GetRasterBand(1).WriteArray(basin.astype(np.uint8))
out = None
