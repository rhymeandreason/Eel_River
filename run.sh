#!/usr/bin/env bash
# 3DEP DEM -> GRASS r.sim.water -> water-depth overlay + terrain for viewer.html.
# Usage: ./run.sh   (uses the `grass` conda env if grass is not already on PATH)
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v grass >/dev/null; then
  exec conda run -n grass --no-capture-output "$0" "$@"
fi

# ---- parameters (override from the environment, e.g. RAIN=25 ./run.sh) ----
export LAT=${LAT:-39.38} LON=${LON:--123.07}  # box centre
export HALF=${HALF:-1000}       # sim box half-width, m  -> 2 x 2 km
export PAD=${PAD:-2000}         # terrain context around the box, m
export RES=${RES:-3}            # sim cell size, m
export TRES=${TRES:-10}         # viewer terrain cell size, m
export RAIN=${RAIN:-50}         # rainfall rate, mm/hr, uniform
export INFIL=${INFIL:-0}        # infiltration, mm/hr (0 = saturated soil)
export MANNING=${MANNING:-0.1}  # Manning's n, one value for the whole box
export MINUTES=${MINUTES:-60}   # storm duration simulated, min
export DMIN=${DMIN:-0.005}      # depths below this (m) are transparent
export NPROCS=${NPROCS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
CRS=EPSG:26910                  # NAD83 / UTM 10N, metres

# ---- stage 2: runs inside a throwaway GRASS project (see bottom) ----
if [[ "${1:-}" == --in-grass ]]; then
  r.in.gdal -o input=data/dem.tif output=dem --overwrite
  g.region raster=dem -p

  # Slope components the SIMWE walkers move along.
  r.slope.aspect elevation=dem dx=dx dy=dy --overwrite

  r.sim.water elevation=dem dx=dx dy=dy \
    rain_value="$RAIN" infil_value="$INFIL" man_value="$MANNING" \
    niterations="$MINUTES" output_step="$MINUTES" random_seed=1 nprocs="$NPROCS" \
    depth=depth discharge=discharge --overwrite

  r.out.gdal input=depth output=data/water_depth_m.tif type=Float32 \
    createopt=COMPRESS=DEFLATE --overwrite
  r.univar -g map=depth > data/depth_stats.txt

  # Colour + alpha: shallow sheet flow faint, channels opaque from 0.3 m.
  r.mapcalc --overwrite "depth_vis = if(depth < $DMIN, null(), depth)"
  r.colors map=depth_vis rules=- <<EOF
$DMIN 199 233 255
0.03 107 174 214
0.1  33 113 181
0.3  8 69 148
1    8 29 88
100  8 29 88
EOF
  r.mapcalc --overwrite <<EOF
w1_r = if(isnull(depth_vis), 0, r#depth_vis)
w2_g = if(isnull(depth_vis), 0, g#depth_vis)
w3_b = if(isnull(depth_vis), 0, b#depth_vis)
w4_a = if(isnull(depth_vis), 0, round(255 * min(1.0, 0.25 + 0.75 * depth_vis / 0.3)))
EOF
  i.group group=water input=w1_r,w2_g,w3_b,w4_a
  r.out.gdal input=water output=data/water_depth_rgba.tif type=Byte \
    createopt=PHOTOMETRIC=RGB,ALPHA=YES,COMPRESS=DEFLATE -c --overwrite
  exit 0
fi

# ---- stage 1: download the DEM ----
mkdir -p data web
read -r CX CY < <(echo "$LON $LAT" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)

fetch_dem() {  # $1 out file, $2 half-width m, $3 cell m
  local x0 y0 n
  x0=$(python -c "print(round($CX/$3)*$3 - $2)")
  y0=$(python -c "print(round($CY/$3)*$3 - $2)")
  n=$(( 2 * $2 / $3 ))
  curl -sSf -o "$1" "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?bbox=$x0,$y0,$((x0 + n * $3)),$((y0 + n * $3))&bboxSR=26910&imageSR=26910&size=$n,$n&format=tiff&pixelType=F32&interpolation=RSP_BilinearInterpolation&f=image"
  gdal_edit -a_srs $CRS "$1"
}
[[ -s data/dem.tif ]]     || fetch_dem data/dem.tif     "$HALF"           "$RES"
[[ -s data/terrain.tif ]] || fetch_dem data/terrain.tif "$((HALF + PAD))" "$TRES"

# ---- stage 2: GRASS ----
grass --tmp-project $CRS --exec "$0" --in-grass

# ---- stage 3: reproject to lon/lat for Cesium ----
gdalwarp -q -overwrite -t_srs EPSG:4326 -r near data/water_depth_rgba.tif data/water_depth_4326.tif
gdal_translate -q -of PNG data/water_depth_4326.tif web/water_depth.png
gdalwarp -q -overwrite -t_srs EPSG:4326 -r bilinear -dstnodata -9999 -of ENVI -ot Float32 \
  data/terrain.tif web/terrain.f32
rm -f web/terrain.f32.aux.xml web/water_depth.png.aux.xml

python - <<'EOF'
import json, os
from osgeo import gdal
gdal.UseExceptions()
def bounds(p):
    d = gdal.Open(p); x0, dx, _, y0, _, dy = d.GetGeoTransform()
    return dict(west=x0, north=y0, east=x0 + dx * d.RasterXSize,
                south=y0 + dy * d.RasterYSize, width=d.RasterXSize, height=d.RasterYSize)
stats = dict(l.split("=") for l in open("data/depth_stats.txt").read().split())
e = os.environ
json.dump(dict(
    terrain=bounds("web/terrain.f32"), water=bounds("data/water_depth_4326.tif"),
    params=dict(lat=float(e["LAT"]), lon=float(e["LON"]), box_m=2 * int(e["HALF"]),
                cell_m=float(e["RES"]), rain_mm_hr=float(e["RAIN"]),
                infil_mm_hr=float(e["INFIL"]), manning_n=float(e["MANNING"]),
                minutes=float(e["MINUTES"]), min_depth_m=float(e["DMIN"])),
    depth_m=dict(max=float(stats["max"]), mean=float(stats["mean"])),
), open("web/meta.json", "w"), indent=1)
EOF
echo "done: web/ has terrain.f32, water_depth.png, meta.json. Serve with: python3 -m http.server -d web 8000"
