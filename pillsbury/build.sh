#!/usr/bin/env bash
# Lake Pillsbury drained: PG&E 2023 bathymetry (from its FERC filing) + USGS 3DEP lidar
# -> one lakebed DEM -> web/ for index.html.   Usage: ./build.sh
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v grass >/dev/null; then
  exec conda run -n grass --no-capture-output "$0" "$@"
fi

export LAT=39.425 LON=-122.945  # lake centre
export HALF=5001                # lidar box half-width, m
export RES=5                    # lakebed DEM cell, m
export WEBRES=8                 # viewer terrain cell, m
CRS=EPSG:26910

# ---- stage 2: GRASS ----
if [[ "${1:-}" == --in-grass ]]; then
  WL_M=$(python -c "import json; print(json.load(open('data/extract.json'))['lidar_water_ft_navd88'] * 1200 / 3937)")
  r.in.gdal -o input=data/lidar_utm.vrt output=lidar --overwrite
  g.region raster=lidar res="$RES" -a
  r.in.gdal -o input=data/footprint_utm_in.tif output=foot --overwrite
  r.null map=foot setnull=0
  r.in.gdal -o input=data/bathy_hue_utm.tif output=bathy_hue --overwrite
  v.import input=data/contours.gpkg output=contours --overwrite

  # Interpolate between the labelled contours, over the surveyed lake only. Slow
  # (~25 min at 5 m), so the result is kept in data/ and reused.
  if [[ -s data/bathy_contour.tif ]]; then
    r.in.gdal -o input=data/bathy_contour.tif output=bathy --overwrite
  else
    g.region vector=contours res="$RES" -a
    v.to.rast input=contours output=crast use=attr attribute_column=navd88_m --overwrite
    r.surf.contour input=crast output=bathy --overwrite
    r.out.gdal -f input=bathy output=data/bathy_contour.tif type=Float32 createopt=COMPRESS=DEFLATE --overwrite
  fi

  # Survey wherever the lidar saw water; lidar everywhere else (it is the measured ground).
  g.region raster=lidar res="$RES" -a
  r.mapcalc --overwrite "lakebed = if(!isnull(foot) && !isnull(bathy) && lidar < $WL_M + 0.1, bathy, lidar)"
  r.mapcalc --overwrite "qa_diff = if(!isnull(foot) && lidar < $WL_M + 0.1, bathy - bathy_hue, null())"
  r.univar -g map=qa_diff > data/qa_contour_vs_ramp.txt
  r.out.gdal -f input=lakebed output=data/lakebed_navd88_m.tif type=Float32 createopt=COMPRESS=DEFLATE --overwrite
  r.mapcalc --overwrite "surveyed = if(isnull(foot), 0, 1)"
  r.out.gdal -f input=surveyed output=data/footprint_utm.tif type=Byte nodata=255 createopt=COMPRESS=DEFLATE --overwrite
  exit 0
fi

mkdir -p data web

# ---- stage 1: inputs ----
# PG&E's 2023 Lake Pillsbury bathymetric survey, FERC eLibrary accession 20240319-5033.
[[ -s data/bathy2023.pdf ]] || curl -sSf -X POST -H "Content-Type: application/json" -o data/bathy2023.pdf \
  -d '{"FileType":"","accession":"","fileid":0,"FileIDAll":"","fileidLst":["2A575A3F-DB71-CDBC-9E8E-8E56B9C00000"],"Islegacy":false}' \
  https://elibrary.ferc.gov/eLibrarywebapi/api/File/DownloadP8File

read -r CX CY < <(echo "$LON $LAT" | gdaltransform -s_srs EPSG:4326 -t_srs $CRS -output_xy)
X0=$(python -c "print(int($CX / 3) * 3 - $HALF)"); Y0=$(python -c "print(int($CY / 3) * 3 - $HALF)")
for i in 0 1; do  # two halves: one 10 km request is over the service's size limit
  f=data/lidar_utm_$i.tif; y=$((Y0 + i * HALF))
  [[ -s $f ]] || { curl -sSf -o $f "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?bbox=$X0,$y,$((X0 + 2 * HALF)),$((y + HALF))&bboxSR=26910&imageSR=26910&size=$((2 * HALF / 3)),$((HALF / 3))&format=tiff&pixelType=F32&interpolation=RSP_BilinearInterpolation&f=image"; gdal_edit -a_srs $CRS $f; }
done
gdalbuildvrt -q -overwrite data/lidar_utm.vrt data/lidar_utm_0.tif data/lidar_utm_1.tif
# the survey's own frame: CA State Plane II, US survey feet
gdalwarp -q -overwrite -t_srs EPSG:2226 -tr 10 10 -r bilinear data/lidar_utm.vrt data/lidar_m_2226.tif
gdal_calc --quiet -A data/lidar_m_2226.tif --calc="A*3937/1200" --NoDataValue=-9999 --outfile=data/lidar_ft.tif --overwrite

python extract.py
# the ramp rasters come out in the sheet's slightly rotated frame; GRASS wants north-up
gdalwarp -q -overwrite -t_srs $CRS -tr 2 2 -r near data/footprint.tif data/footprint_utm_in.tif
gdalwarp -q -overwrite -t_srs $CRS -tr 2 2 -r bilinear data/bathy_hue.tif data/bathy_hue_utm.tif

# ---- stage 2 ----
grass --tmp-project $CRS --exec "$0" --in-grass

# ---- stage 3: web ----
python basin.py
gdalwarp -q -overwrite -t_srs EPSG:4326 -tr $(python -c "print($WEBRES/111320, $WEBRES/111320)") -r bilinear \
  -dstnodata -9999 -of ENVI -ot Float32 data/lakebed_navd88_m.tif web/lakebed.f32
gdalwarp -q -overwrite -t_srs EPSG:4326 -te $(python - <<'EOF'
from osgeo import gdal
d = gdal.Open("web/lakebed.f32"); x0, dx, _, y0, _, dy = d.GetGeoTransform()
print(x0, y0 + dy * d.RasterYSize, x0 + dx * d.RasterXSize, y0, "-ts", d.RasterXSize, d.RasterYSize)
EOF
) -r near -of ENVI -ot Byte data/basin.tif web/basin.u8
rm -f web/*.aux.xml

python - <<'EOF'
import json
from osgeo import gdal
gdal.UseExceptions()
d = gdal.Open("web/lakebed.f32"); x0, dx, _, y0, _, dy = d.GetGeoTransform()
ex = json.load(open("data/extract.json"))
qa = dict(l.split("=") for l in open("data/qa_contour_vs_ramp.txt").read().split())
json.dump(dict(
    grid=dict(west=x0, north=y0, east=x0 + dx * d.RasterXSize, south=y0 + dy * d.RasterYSize,
              width=d.RasterXSize, height=d.RasterYSize),
    pge_minus_navd88_ft=ex["pge_minus_navd88_ft"],
    full_pool_pge_ft=1910, spill_crest_pge_ft=1900, survey_water_pge_ft=1897.3,
    lidar_water_pge_ft=ex["lidar_water_ft_navd88"] + ex["pge_minus_navd88_ft"],
    contour_vs_ramp_m=dict(mean=float(qa["mean"]), sd=float(qa["stddev"])),
    ramp_lidar_agree=ex["ramp_lidar_agree"],
), open("web/meta.json", "w"), indent=1)
EOF
echo "done: serve with  python3 -m http.server -d web 8001"
