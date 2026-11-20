# NWS Hazard Merge

Merges the NWS national county/zone hazard file with the NDFD gridded hazard
grids (days 1–3 and 4–7). It writes one shapefile and one GeoJSON, sorted for
hazard-priority drawing.

## Inputs (downloaded on every run)

| Source | URL |
|---|---|
| National WWA (county/zone) | `https://tgftp.nws.noaa.gov/SL.us008001/DF.sha/DC.cap/DS.WWA/current_all.tar.gz` |
| NDFD hazards, days 1–3 (hourly) | `.../DC.ndfd/AR.conus/VP.001-003/ds.wwa.bin` |
| NDFD hazards, days 4–7 (6-hourly) | `.../DC.ndfd/AR.conus/VP.004-007/ds.wwa.bin` |
| Priority and colors | `config/hazard_priority.csv`, built from https://www.weather.gov/help-map (March 10, 2025 list) |

## How the merge works

1. **National features** keep all original attributes. Their valid window is
   `VALID_STRT` = ONSET (or ISSUANCE if missing) and `VALID_END` = ENDS (or
   EXPIRATION if missing). Features with a VTEC action of CAN, EXP or UPG are
   dropped, and so are already-expired ones (`--keep-expired` turns this off).
2. **Only CONUS features are checked against NDFD** (feature entirely inside
   `--conus-bbox`). Alaska, Hawaii, Puerto Rico, Guam and others pass through
   with their original zone geometry, unchanged.
3. **Matching.** The NDFD grid stores strings such as `CF.W^FA.A` per cell. A
   national feature matches when its hazard code (PHENOM.SIG) is in NDFD cells
   inside the zone during any NDFD time step that overlaps its valid window
   (±`--time-tolerance-min`, default 60). The result is shown in `GEOM_SRC`:
   * `ZONE`: no NDFD match. Zone geometry is kept.
   * `ZONE_FULL`: NDFD covers at least 95% of the zone (`--full-zone-frac`).
     The exact zone outline is kept.
   * `NDFD_CLIP`: NDFD covers only part of the zone. The geometry becomes the
     zone outline minus the zone cells with no hazard, so zone edges and
     coastlines stay crisp and NDFD supplies the interior detail.
     `NDFD_FRAC` is the fraction of zone cells covered.
4. **NDFD-only hazards.** Grid cells and times that no national feature of the
   same code explains are turned into polygons (`GEOM_SRC = NDFD_ONLY`,
   `SOURCE = NDFD`). Cells with the same on/off timeline are grouped, so each
   feature has a single continuous `VALID_STRT`–`VALID_END`. Features with the
   same hazard never overlap in time at the same place. These are often
   offshore or outer marine zones and days 4–7 watches that aren't issued yet.
5. **Slivers.** Polygon parts and holes smaller than `--min-part-cells` (2
   cells, about 13 km²) are dropped. NDFD-derived edges are simplified by
   `--simplify-m` (500 m, 0 = off).
6. **Draw order.** Features are written from lowest priority (largest number)
   to highest priority (1). Renderers that draw in file order (QGIS, ArcGIS,
   Leaflet, Mapbox, …) therefore put the highest priority on top. `DRAW_ORDR`
   holds the same order (1 = drawn first/bottom). If you use symbol levels, sort
   by `PRIORITY` descending.

NDFD time-step convention: each grid is treated as valid from its valid time
until the next grid's valid time. The last 6-hourly grid lasts 6 hours.

## Output fields

`HAZ_CODE, HAZ_NAME, PRIORITY, COLOR (hex), DRAW_ORDR, VALID_STRT, VALID_END (UTC ISO-8601), GEOM_SRC, NDFD_FRAC, SOURCE`,
followed by the original national fields
`CAP_ID, MSG_TYPE, PROD_TYPE, VTEC, PHENOM, SIG, WFO, EVENT, ISSUANCE, EXPIRATION, ONSET, ENDS, URL`
(left empty for NDFD-only features). CRS is EPSG:4326. Files are replaced
atomically at the end of each run.

## Install (Linux)

GDAL comes bundled in the `rasterio`/`pyogrio` wheels, so no system GDAL is
needed. Requires Python 3.10 or newer.

```bash
cd nws_hazard_merge
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run

```bash
.venv/bin/python hazard_merge.py                    # download + merge -> output/merged_hazards.*
.venv/bin/python hazard_merge.py --no-download      # reuse files in work/
.venv/bin/python hazard_merge.py --help             # all tuning options
```

## Cron

`run_hazard_merge.sh` takes a lock (overlapping runs are skipped) and logs to
`logs/hazard_merge.log`. A run takes about 30 seconds. Example, every 10 minutes:

```
*/10 * * * * /opt/nws_hazard_merge/run_hazard_merge.sh --outdir /data/hazards
```

## Maintenance

* **Priority list changes:** edit `config/hazard_priority.csv` (priority,
  hazard name, space-separated VTEC codes, color). Hazards not in the table get
  priority 999 and gray, and the log warns about them.
* **GRIB decoding:** the grid is read with GDAL. GDAL handles the NDFD
  alternating-row scan order; eccodes does not. The hazard string table is
  decoded from the GRIB2 Local Use section by `_decode_wx_table`.
