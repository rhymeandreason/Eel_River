"""Check the lakebed DEM against PG&E's own 2023 area-capacity table (report 026.11-24.2,
in the same FERC filing). Below 1897 ft the table comes from the same sonar survey, so a
DEM built correctly from its contours should reproduce it."""
import json
import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()
# PG&E ft: (acres, acre-ft). 2023 survey up to 1896.7; above that the 2023 survey missed
# the north end, so those rows are the 2016 survey's.
TABLE = {1851.7: (390, 5768), 1861.7: (573, 10521), 1871.7: (804, 17313),
         1881.7: (1153, 27057), 1891.7: (1477, 40365), 1896.7: (1556, 47959),
         1900.0: (1875, 56322), 1910.0: (2275, 76876)}

off = json.load(open("data/extract.json"))["pge_minus_navd88_ft"]
ds = gdal.Open("data/lakebed_navd88_m.tif")
z = ds.ReadAsArray().astype(float)
nd = ds.GetRasterBand(1).GetNoDataValue()
z[z == nd] = np.nan
cell = ds.GetGeoTransform()[1] ** 2
foot = gdal.Open("data/footprint_utm.tif").ReadAsArray() == 1
basin = gdal.Open("data/basin.tif").ReadAsArray() == 1
seed = np.unravel_index(np.nanargmin(np.where(foot, z, np.nan)), z.shape)
print(f"deepest surveyed point {(z[seed] / 0.3048006 + off):.1f} ft PG&E (table's zero capacity: 1817.7)")
print(f"{'level':>7} {'acres':>14} {'acre-ft':>18}")
for lvl, (ac, af) in TABLE.items():
    m = (lvl - off) * 0.3048006
    lab, _ = ndimage.label((z < m) & basin)
    wet = lab == lab[seed]
    a = wet.sum() * cell / 4046.86
    v = np.nansum((m - z)[wet]) * cell / 1233.48
    print(f"{lvl:7.1f} {a:6.0f} vs {ac:5d} {v:8.0f} vs {af:6d}  ({100 * (v / af - 1):+.1f}%)")
