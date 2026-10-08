# Eel River runoff: GRASS r.sim.water → CesiumJS

See also [`pillsbury/`](pillsbury/README.md): Lake Pillsbury's floor from PG&E's 2023 survey, drained in a 3D viewer.

A 2 × 2 km box of the upper Eel River, about 4 km upstream (east) of Cape Horn Dam, Mendocino County, CA. One uniform storm falls on a USGS lidar DEM, GRASS routes the water over the surface, and a web page drapes the resulting water depth on 3D terrain.

```
run.sh          download → GRASS → export (the whole pipeline)
web/index.html  the viewer
web/            what the viewer loads (meta.json, terrain.f32, water_depth.png)
data/           intermediate GeoTIFFs, kept for QGIS etc.
```

## Run it

```bash
./run.sh
```

```bash
python3 -m http.server -d web 8000
```

Then open http://localhost:8000. The page fetches files, so it will not work from `file://`.

`run.sh` finds GRASS in the `grass` conda env (`conda create -n grass -c conda-forge grass gdal`, ~670 MB). On an M-series Mac the sim takes about 7 minutes. Every parameter is an environment variable at the top of the script, e.g. `RAIN=25 MINUTES=30 ./run.sh`. The DEM is only downloaded once; delete `data/*.tif` to fetch again (needed if you change `LAT`, `LON`, `HALF` or `RES`).

## What each step does

### 1. Download the DEM (`run.sh` stage 1)

The USGS 3DEP elevation service (`elevation.nationalmap.gov/.../3DEPElevation/ImageServer/exportImage`) serves the best bare-earth elevation it has for any box, resampled on request. No key needed. Here it is lidar-derived (roads and gullies are resolved). Fetched twice:

| File | Extent | Cell | Used for |
| --- | --- | --- | --- |
| `data/dem.tif` | 2 × 2 km | 3 m | the simulation |
| `data/terrain.tif` | 6 × 6 km | 10 m | the viewer's 3D terrain, with 2 km of context around the box |

Both are requested in NAD83 / UTM zone 10N (EPSG:26910) so cells are square metres, which `r.slope.aspect` and `r.sim.water` assume. The box centre (39.38, -123.07) is converted to UTM with `gdaltransform` and snapped to the cell grid.

### 2. GRASS (`run.sh` stage 2)

`grass --tmp-project EPSG:26910 --exec run.sh --in-grass` makes a throwaway GRASS project, reruns the same script inside it, and deletes the project afterwards. Inside:

1. **`r.in.gdal`** imports the DEM; `g.region raster=dem` sets the computational region to its extent and 3 m cells.
2. **`r.slope.aspect dx=dx dy=dy`** writes the slope's east-west and north-south components (m/m). These are the direction and steepness water moves along.
3. **`r.sim.water`** (SIMWE) solves shallow overland flow by path sampling: it releases "walkers", each a small parcel of water, which move downslope by the dx/dy gradient plus a diffusion term, at a speed set by Manning's equation. Counting walkers per cell gives depth. Settings:

   | Parameter | Value | Meaning |
   | --- | --- | --- |
   | `rain_value` | 50 mm/hr | uniform over the box. An intense storm for this region |
   | `infil_value` | 0 mm/hr | saturated soil, so all rain runs off. Worst case |
   | `man_value` | 0.1 | one Manning's n for everything |
   | `niterations` | 60 min | storm length; depth is reported at the end |
   | `random_seed` | 1 | same walkers every run, so results are reproducible |

   Outputs: `depth` (m) and `discharge` (m³/s).
4. **`r.out.gdal`** writes the raw depth as `data/water_depth_m.tif` (Float32, metres). This is the file to analyse.
5. **Colourize.** Depths under 5 mm become null. `r.colors` sets a blue ramp from 5 mm to 1 m. `r.mapcalc` splits the colour table into red, green and blue bands (`r#`, `g#`, `b#`) and builds an alpha band that runs from 25% for a film of sheet flow up to opaque at 30 cm. `i.group` + `r.out.gdal` write `data/water_depth_rgba.tif`, an RGBA GeoTIFF with transparent dry cells.

### 3. Export for the web (`run.sh` stage 3)

Cesium drapes an image on a lon/lat rectangle, so both layers are reprojected to WGS84 (EPSG:4326) with `gdalwarp`:

- `web/water_depth.png`: the RGBA overlay (nearest-neighbour, so colours are not blended).
- `web/terrain.f32` + `.hdr`: the 10 m context DEM as a raw little-endian Float32 grid (GDAL's ENVI format).
- `web/meta.json`: both rectangles, the sim parameters, and max/mean depth. The viewer reads its numbers from here, never types them.

### 4. Viewer (`web/index.html`)

CesiumJS 1.146 from jsDelivr, no Cesium ion account:

- **Terrain**: a `CustomHeightmapTerrainProvider` samples `terrain.f32` (bilinear) for every tile Cesium asks for, so the water is draped on the same lidar surface it was simulated on. Outside the 6 km context the ground sits flat at the DEM's lowest height.
- **Imagery**: Esri World Imagery tiles (keyless).
- **Water**: `water_depth.png` via `SingleTileImageryProvider` over its rectangle. Being imagery, it drapes on the terrain and cannot float or sink.
- **Controls**: an opacity slider (sets the layer's `alpha`) and an optional 2× vertical exaggeration.

## Reading the result, and what it leaves out

- **The main Eel channel only carries rain that fell inside the box.** There is no inflow from upstream (Lake Pillsbury and above) and no baseflow, so its depth here is the local storm's contribution, far below a real flood stage. Adding inflow means a `rain` raster with a source cell at the upstream edge, or a 1D/2D hydraulic model.
- **The river bed is the water surface.** Lidar does not see through water, so the "bed" in the DEM is the low-flow water surface on the day of the flight.
- **The deepest cells (max ~10.7 m) are ponding behind road fills.** Lidar sees the road embankment but not the culvert under it, so the sim fills the hollow upstream like a dam. Only a few dozen cells exceed 3 m. Hydro-conditioning the DEM (burning culverts through roads) would drain them.
- **One Manning's n.** Real values span ~0.03 (open channel) to 0.4+ (forest litter). Using one value makes hillslope flow too fast and channel flow too slow. `r.sim.water man=<raster>` takes a map, e.g. from NLCD land cover.
- **Heights are NAVD88, not the ellipsoid Cesium expects.** Everything shown is from the same DEM, so it is consistent; it would sit ~30 m high next to Cesium World Terrain (geoid is ~-32 m here).
