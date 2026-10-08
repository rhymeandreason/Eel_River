#!/usr/bin/env bash
# The rivers across the drained lakebed. Run after build.sh.   Usage: ./rivers.sh
#   1. take the dam out of the DEM (breach.py)
#   2. route water over it in GRASS, with the drainage area from outside the 10 km box
#      carried in from a 30 m DEM of the whole upper Eel
#   3. turn three measured flows into water depths (rivers.py)
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v grass >/dev/null; then
  exec conda run -n grass --no-capture-output "$0" "$@"
fi
CRS=EPSG:26910
export STREAM_KM2=${STREAM_KM2:-1}   # drainage area at which a channel starts

if [[ "${1:-}" == --in-grass ]]; then
  read -r DX DY < <(python -c "import json; d = json.load(open('data/dam.json')); print(d['utm_x'], d['utm_y'])")

  # whole upper Eel at 30 m: how much land drains to the dam, and from where
  r.in.gdal -o input=data/region30.tif output=region30 --overwrite
  g.region raster=region30
  r.watershed -s elevation=region30 accumulation=acc30 drainage=dir30 --overwrite
  r.mapcalc --overwrite "near = if(sqrt((x() - $DX)^2 + (y() - $DY)^2) < 600, abs(acc30), null())"
  read -r OX OY < <(r.stats -gn near | sort -k3 -g | tail -1 | cut -d' ' -f1,2)
  r.water.outlet input=dir30 output=basin30 coordinates="$OX,$OY" --overwrite
  r.stats -an basin30 | awk '{print "upper_eel_km2=" $2 / 1e6}' > data/rivers_grass.txt
  r.mapcalc --overwrite "edge = if(!isnull(basin30) && (row() < 3 || row() > nrows() - 3 || col() < 3 || col() > ncols() - 3), 1, null())"
  [[ -z "$(r.stats -n edge)" ]] || { echo "ERROR: the upper Eel runs off region30.tif; widen its box" >&2; exit 1; }

  # Route over the whole upper Eel at 10 m with the drained lakebed patched in, so every
  # river reaches the lake on its own (a box around the lake cuts rivers at its edges, and
  # r.watershed treats edge cells as outlets).
  r.in.gdal -o input=data/region10.vrt output=region10 --overwrite
  r.in.gdal -o input=data/drained_navd88_m.tif output=drained --overwrite
  g.region raster=region10
  r.resamp.stats input=drained output=drained10 method=average --overwrite
  r.patch input=drained10,region10 output=dem10 --overwrite
  r.watershed -s elevation=dem10 accumulation=acc10 drainage=dir10 --overwrite
  r.stream.extract elevation=dem10 accumulation=acc10 threshold=$((STREAM_KM2 * 10000)) \
    stream_raster=streams --overwrite
  # rivers.py needs only the lake's box
  g.region raster=drained res=10 -a
  for m in dem10 acc10 dir10 streams; do
    r.out.gdal -f input=$m output=data/$m.tif type=Float64 createopt=COMPRESS=DEFLATE --overwrite
  done
  exit 0
fi

python breach.py

# 30 m DEM of the upper Eel (the basin above Scott Dam reaches ~40 km north-east)
[[ -s data/region30.tif ]] || {
  read -r X0 Y0 < <(echo "-123.10 39.10" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)
  read -r X1 Y1 < <(echo "-122.60 39.70" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)
  X0=${X0%.*} Y0=${Y0%.*} X1=${X1%.*} Y1=${Y1%.*}
  curl -sSf -o data/region30.tif "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?bbox=$X0,$Y0,$X1,$Y1&bboxSR=26910&imageSR=26910&size=$(((X1 - X0) / 30)),$(((Y1 - Y0) / 30))&format=tiff&pixelType=F32&interpolation=RSP_BilinearInterpolation&f=image"
  gdal_edit -a_srs $CRS data/region30.tif
}

# 10 m DEM over the upper Eel, in tiles under the service's size limit
if [[ ! -s data/region10.vrt ]]; then
  read -r X0 Y0 < <(echo "-123.02 39.18" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)
  read -r X1 Y1 < <(echo "-122.66 39.66" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)
  X0=$(( ${X0%.*} / 10 * 10 )) Y0=$(( ${Y0%.*} / 10 * 10 )) X1=${X1%.*} Y1=${Y1%.*}
  T=20000
  for ((x = X0; x < X1; x += T)); do for ((y = Y0; y < Y1; y += T)); do
    f=data/r10_${x}_${y}.tif
    [[ -s $f ]] || { curl -sSf -o $f "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?bbox=$x,$y,$((x + T)),$((y + T))&bboxSR=26910&imageSR=26910&size=$((T / 10)),$((T / 10))&format=tiff&pixelType=F32&interpolation=RSP_BilinearInterpolation&f=image"; gdal_edit -a_srs $CRS $f; }
  done; done
  gdalbuildvrt -q data/region10.vrt data/r10_*.tif
fi

grass --tmp-project $CRS --exec "$0" --in-grass
python rivers.py
