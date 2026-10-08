# Removing Scott Dam

Since 1922, Scott Dam has blocked the upper Eel River, holding back Lake Pillsbury and cutting salmon and steelhead off from the headwaters. PG&E now plans to surrender its Potter Valley Project and take the dam down. Federal and state environmental reviews are under way.

This is a 3D viewer of what removal would expose. The ground is the floor of Lake Pillsbury as PG&E surveyed it in July 2023, joined to USGS lidar. A water-level slider runs from full pool (1910 ft) down to the lake floor. With the lake gone, the rivers that feed it still cross the old bed, shown at three measured flows.

```
build.sh       download → extract contours from the PDF → GRASS surface → merge → web/
extract.py     the PDF part: georeference, datum, and an elevation for every contour
basin.py       where water can stand (cuts the canyon below the dam off)
qa_capacity.py checks the result against PG&E's own area-capacity table
rivers.sh      the rivers: cut the dam out (breach.py), route the upper Eel, depths (rivers.py)
web/index.html the viewer
```

```bash
./build.sh
```

```bash
./rivers.sh
```

```bash
python3 -m http.server -d web 8001
```

Elevations are in **PG&E datum feet** (the datum every PG&E and FERC document for this dam uses). PG&E datum = NAVD88 + 78.55 ft here.

## Where the lake floor comes from

There is no published grid of the lake floor. PG&E's 2023 bathymetric survey report is public on FERC eLibrary (accession 20240319-5033), and its Sheet 1 is an AutoCAD export: 5 ft contours as vector paths over a colour-ramped depth image, in CA State Plane II feet. `extract.py` turns that into labelled contour lines:

1. **Contours.** About 15,000 path fragments on the PDF layer "Surface 2023", chained into 10,045 lines. The darker lines are the 25 ft majors.
2. **Position.** Scale from the sheet (1 pt = 27.78 ft, from the viewport and the scale bar). The offset comes from sliding the survey's outline over the lake in the 3DEP lidar. Shift, rotation and scale are then refined so the contours above the lidar's water line run as level as possible on the lidar. The fitted rotation, −0.59°, matches the State Plane grid convergence here (about 0.6°), so the sheet appears to be drawn to true north.
3. **Datum.** Contours on dry lidar ground must fall on 5 ft multiples, and the majors on 25 ft multiples. That gives PG&E = NAVD88 + 78.55 ft. Independent check: the report's own capacity table puts PG&E datum 81.7 ft above the USGS (NGVD29) datum, and NGVD29 → NAVD88 is about +3 ft in this area.
4. **Elevations.** The contours have no printed labels. Dry contours take their elevation from the lidar. Submerged ones take it from the colour of the depth image under the line: contour hues fall in tight clusters, with four minor clusters between each pair of majors, which pins each cluster to a 5 ft level (`HUE_LADDER`). Where both readings exist they agree on 14 of 15 lines.

`build.sh` then:

5. **Interpolates** the contours with GRASS `r.surf.contour` at 5 m (UTM 10N, NAVD88 m). It is slow (about 25 min), so the result is cached in `data/bathy_contour.tif`.
6. **Merges.** The survey surface is used wherever the lidar saw water; the lidar is used everywhere else. The lidar was flown with the lake drawn down to 1888 ft, so it already measures the upper lakebed, including the delta flats at the north end that the 2023 boats could not reach.
7. **Draws the basin** (`basin.py`). The lidar is bare earth, so Scott Dam is not in it, and water above ~1891 ft would run straight down the Eel. The script finds the canyon below the dam (ground deeper than the lake floor that the lake does not reach), plugs a 150 m circle where it comes closest to the lake, and keeps the cells below 1910 ft that connect to the lake bottom. The plug lands at 39.407, −122.959, the dam site.

## The rivers on the drained bed (`rivers.sh`)

8. **Take the dam out** (`breach.py`). The lidar still has a ridge about 1890 ft high at the dam site, so a drained basin would hold a 75 ft pond. The script cuts an 860 m notch on a straight grade from the lake floor (1815 ft) to the riverbed below the dam (1792 ft). It only lowers the ground. The highest ground cut is 1899 ft, which is the dam ridge, not a hillside. The real removal design may excavate differently.
9. **Route the water.** GRASS `r.watershed` runs over a 10 m DEM of the whole upper Eel, with the drained lakebed patched in, so every river reaches the lake on its own. A 30 m run delineates the basin above the dam at 289 sq mi, against USGS's 290 for gauge 11470500. At 10 m, 282 sq mi of it drains out through the notch. Channels start at 1 km² of drainage (`r.stream.extract`).
10. **Pick the flows.** The gauge below Scott Dam measures released water, not the river, so the flows come from USGS 11473900, Middle Fork Eel near Dos Rios (745 sq mi, no dam, 1965–2025), scaled by drainage area. Its median annual flood is 46 cfs per sq mi, against 41 at the Scott Dam gauge. At the dam site that gives:

    | Flow | From | cfs |
    | --- | --- | --- |
    | Late summer | median August day | 7 |
    | Typical winter | median January–February day | 756 |
    | 2-year flood | median annual peak | 12,938 |

11. **Depths** (`rivers.py`). HAND with a synthetic rating curve per reach, the method behind NOAA's National Water Model flood maps. HAND is each cell's height above the channel cell it drains to. Filling a reach's cells to stage h gives a flow area and wetted perimeter, so Manning's equation (n = 0.05) gives Q(h). Inverting at the reach's flow gives its stage, and depth = stage − HAND.

What to read into it:

- **It is the first winter, not the river that settles.** Water follows the 2023 sediment surface. The flood map shows where water would spread across an un-incised bed. Within a few seasons the rivers cut channels into the fill and much of the 21 million cubic yards moves downstream. No sediment is moved here.
- **Ponds are pits in the 2023 surface**, filled to the river's stage. Some are real hollows; some are interpolation between contours.
- **Stage is set reach by reach**, so it can step where two reaches meet. That is why the winter water sometimes reads as pools rather than one ribbon.
- **10 m cells.** The late-summer creek is narrower than one cell, so its map shows where it runs, not how wide it is.

## How good is it

`qa_capacity.py` floods the DEM and compares the result with the area-capacity table in the same FERC filing. Rows up to 1896.7 ft are PG&E's 2023 survey. Above that the 2023 survey missed the north end, so those rows are PG&E's 2016 survey.

| Level (ft) | Area, ours vs PG&E (acres) | Volume, ours vs PG&E (acre-ft) |
| --- | --- | --- |
| 1851.7 | 389 vs 390 | 5,272 vs 5,768 (−8.6%) |
| 1871.7 | 806 vs 804 | 16,652 vs 17,313 (−3.8%) |
| 1891.7 | 1,565 vs 1,477 | 38,857 vs 40,365 (−3.7%) |
| 1910.0 | 2,271 vs 2,275 | 73,802 vs 76,876 (−4.0%) |

Areas match closely. Volume runs 2–9% low, mostly near the bottom: the lowest contour is 1815 ft, so pockets below it (the table's zero is 1817.7) are filled flat, and `r.surf.contour` rounds off the channel floors.

## What it does and doesn't show

- **The drained surface is the sediment, not the old valley.** About 21 million cubic yards have settled behind Scott Dam since 1921, burying the 1921 valley floor. After removal the Eel cuts down into this surface, and much of the fill moves downstream over years. Showing that needs the pre-dam (1922) surface, which only McBain Associates has digitized, and a sediment-transport model.
- **The deep end is the least certain part.** The colour ramp is compressed in the blues (about 4° of hue per 5 ft below 1850 ft), so labels there could be off by one contour. The capacity table is the check that they are not off by much.
- **The bare-ground colour is a hillshade of the DEM, not a photo.** The satellite imagery under it shows the lake full.
- **The dam is shown removed.** The viewer's terrain has the notch cut through the dam site. At high lake levels the plug stops the water there, but no dam is drawn.
- **Two dates are mixed.** The deep floor is from July 2023; the upper lakebed and north flats are from whenever 3DEP flew this area. PG&E's report measured 1–9 ft of new deposition in the arms between 2016 and 2023, so the flats have kept changing.
