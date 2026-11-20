#!/usr/bin/env python3
"""
Merge NWS county/zone hazards (national WWA shapefile) with NDFD gridded hazards.

Matching rules
--------------
* A national (county/zone) hazard matches NDFD when the NDFD grid contains the
  same VTEC hazard code (PHENOM.SIG, e.g. "WS.W") at one or more NDFD time steps
  that overlap the hazard's valid window [VALID_STRT, VALID_END], inside the
  county/zone polygon.
* Matched features keep every original attribute; only the geometry changes, to
  (zone polygon) INTERSECT (NDFD footprint of that code during the valid window).
* NDFD hazard areas that no national feature explains (same code, same time) are
  converted from grid cells to polygons and added as new features. Their
  VALID_STRT / VALID_END come from the NDFD time steps.
* Output is sorted so the highest priority hazard (lowest PRIORITY number) is
  written LAST, so it draws on top in GIS and web-map renderers that draw
  features in file order. DRAW_ORDR gives the same order explicitly.

Run `python hazard_merge.py --help` for options.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import rasterio.windows
from shapely.geometry import MultiPolygon, Polygon, box, shape
from shapely.ops import unary_union

LOG = logging.getLogger("hazard_merge")

WWA_URL = "https://tgftp.nws.noaa.gov/SL.us008001/DF.sha/DC.cap/DS.WWA/current_all.tar.gz"
NDFD_URLS = {
    "001-003": "https://tgftp.nws.noaa.gov/SL.us008001/ST.opnl/DF.gr2/DC.ndfd/AR.conus/VP.001-003/ds.wwa.bin",
    "004-007": "https://tgftp.nws.noaa.gov/SL.us008001/ST.opnl/DF.gr2/DC.ndfd/AR.conus/VP.004-007/ds.wwa.bin",
}
HERE = Path(__file__).resolve().parent
DEFAULT_PRIORITY_CSV = HERE / "config" / "hazard_priority.csv"
UNKNOWN_PRIORITY = 999
UNKNOWN_COLOR = "#C0C0C0"
NDFD_NODATA = 9999


# --------------------------------------------------------------------------- #
# Download helpers
# --------------------------------------------------------------------------- #
def download(url: str, dest: Path, retries: int = 4, timeout: int = 120) -> Path:
    """Download url to dest atomically, retrying with backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "nws-hazard-merge/1.0"})
    for attempt in range(1, retries + 1):
        try:
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f)
            os.replace(tmp, dest)
            LOG.info("downloaded %s (%d bytes)", url, dest.stat().st_size)
            return dest
        except Exception as exc:  # noqa: BLE001
            if attempt == retries:
                raise
            wait = 2 ** attempt
            LOG.warning("download failed (%s), retry %d in %ds", exc, attempt, wait)
            time.sleep(wait)
    return dest


def extract_shapefile(tgz: Path, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tgz) as tf:
        try:
            tf.extractall(outdir, filter="data")
        except TypeError:  # Python < 3.12
            for m in tf.getmembers():
                if m.isfile() and not os.path.isabs(m.name) and ".." not in Path(m.name).parts:
                    tf.extract(m, outdir)
    shps = sorted(outdir.rglob("*.shp"))
    if not shps:
        raise RuntimeError(f"no shapefile found in {tgz}")
    return shps[0]


# --------------------------------------------------------------------------- #
# Priority / color table
# --------------------------------------------------------------------------- #
@dataclass
class HazardInfo:
    name: str
    priority: int
    color: str


def load_priority_table(path: Path) -> tuple[dict[str, HazardInfo], dict[str, HazardInfo]]:
    by_code, by_name = {}, {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            info = HazardInfo(row["hazard"].strip(), int(row["priority"]), row["hex"].strip())
            by_name[info.name.lower()] = info
            for code in row["vtec_codes"].split():
                by_code.setdefault(code.upper(), info)
    return by_code, by_name


def lookup_hazard(code, prod_type, by_code, by_name) -> HazardInfo:
    if isinstance(code, str) and code.upper() in by_code:
        return by_code[code.upper()]
    if isinstance(prod_type, str) and prod_type.strip().lower() in by_name:
        return by_name[prod_type.strip().lower()]
    name = prod_type if isinstance(prod_type, str) and prod_type else (code or "Unknown")
    return HazardInfo(name, UNKNOWN_PRIORITY, UNKNOWN_COLOR)


# --------------------------------------------------------------------------- #
# NDFD GRIB2 hazard grid
# --------------------------------------------------------------------------- #
def _grib_messages(path: Path):
    """Yield (valid_time_utc, [hazard strings]) for each GRIB2 message.

    The NDFD hazard (WWA) grid stores an integer index per cell. The index refers
    to a table of strings such as "CF.W^FA.A" packed into the GRIB2 Local Use
    section (Section 2): template byte, 2-byte group count, 4-byte number of
    characters, 4-byte reference value, 2 scale bytes, bits-per-character,
    type byte, then the bit-packed NUL-separated characters.
    """
    data = path.read_bytes()
    pos = 0
    while True:
        pos = data.find(b"GRIB", pos)
        if pos < 0:
            return
        total = int.from_bytes(data[pos + 8:pos + 16], "big")
        msg = data[pos:pos + total]
        off, table, valid = 16, None, None
        ref = None
        while off < len(msg) - 4:
            if msg[off:off + 4] == b"7777":
                break
            length = int.from_bytes(msg[off:off + 4], "big")
            sec = msg[off + 4]
            body = msg[off + 5:off + length]
            if sec == 1:
                y = int.from_bytes(body[7:9], "big")
                ref = datetime(y, body[9], body[10], body[11], body[12], body[13], tzinfo=timezone.utc)
            elif sec == 2:
                table = _decode_wx_table(body)
            elif sec == 4:
                tmpl = int.from_bytes(body[2:4], "big")
                if tmpl != 0:
                    raise RuntimeError(f"unexpected product template {tmpl}")
                unit = body[12]
                fct = int.from_bytes(body[13:17], "big")
                mult = {0: 60, 1: 3600, 2: 86400, 10: 3 * 3600, 11: 6 * 3600, 12: 12 * 3600, 13: 1}[unit]
                valid = ref + timedelta(seconds=fct * mult)
            off += length
        if table is None or valid is None:
            raise RuntimeError(f"could not decode hazard table/time in {path}")
        yield valid, table
        pos += total


def _decode_wx_table(body: bytes) -> list[str]:
    ngroups = int.from_bytes(body[1:3], "big")
    if body[0] != 1 or ngroups != 1:
        raise RuntimeError(f"unsupported NDFD local-use section (template {body[0]}, groups {ngroups})")
    nchars = int.from_bytes(body[3:7], "big")
    refval = int(np.frombuffer(body[7:11], ">f4")[0])
    nbits = body[13]
    bits = np.unpackbits(np.frombuffer(body[15:], np.uint8))[: nchars * nbits].reshape(nchars, nbits)
    vals = bits.dot(1 << np.arange(nbits - 1, -1, -1)) + refval
    text = bytes(vals.astype(np.uint8)).decode("ascii", "replace")
    return [s for s in text.split("\0")]


@dataclass
class NDFDGrid:
    times: list[datetime]          # step start (valid time)
    ends: list[datetime]           # step end (next valid time)
    index: np.ndarray              # (T, ny, nx) uint8/uint16 table index (nodata -> 0)
    tables: list[list[set[str]]]   # per step: per index -> set of hazard codes
    transform: rasterio.Affine
    crs: rasterio.crs.CRS
    codes: set[str]


def read_ndfd(paths: list[Path], last_step_hours: float = 6.0) -> NDFDGrid:
    steps = {}  # valid time -> (path, band, table)
    transform = crs = shape_ = None
    for p in paths:
        with rasterio.open(p) as ds:
            if transform is None:
                transform, crs, shape_ = ds.transform, ds.crs, ds.shape
            elif ds.transform != transform or ds.shape != shape_:
                raise RuntimeError(f"{p} has a different grid than the first NDFD file")
            for band, (valid, table) in enumerate(_grib_messages(p), start=1):
                if band > ds.count:
                    raise RuntimeError(f"{p}: GRIB message count mismatch")
                if valid not in steps:  # earlier file (days 1-3) wins on overlap
                    steps[valid] = (p, band, table)
    times = sorted(steps)
    ends = times[1:] + [times[-1] + timedelta(hours=last_step_hours)]
    maxlen = max(len(steps[t][2]) for t in times)
    dtype = np.uint8 if maxlen < 255 else np.uint16
    index = np.zeros((len(times),) + shape_, dtype=dtype)
    tables, codes = [], set()
    open_ds = {}
    try:
        for i, t in enumerate(times):
            p, band, table = steps[t]
            ds = open_ds.get(p) or open_ds.setdefault(p, rasterio.open(p))
            a = ds.read(band)
            a = np.where((a == NDFD_NODATA) | ~np.isfinite(a), 0, a)
            index[i] = a.astype(dtype)
            parsed = []
            for s in table:
                s = s.strip()
                hz = set()
                if s and not s.startswith("<"):
                    for h in s.split("^"):
                        h = h.split(":")[0].strip().upper()
                        if len(h) == 4 and h[2] == ".":
                            hz.add(h)
                parsed.append(hz)
            tables.append(parsed)
            codes |= set().union(*parsed) if parsed else set()
    finally:
        for ds in open_ds.values():
            ds.close()
    LOG.info("NDFD: %d steps %s -> %s, hazards: %s", len(times), times[0], ends[-1], sorted(codes))
    return NDFDGrid(times, ends, index, tables, transform, crs, codes)


def code_mask(grid: NDFDGrid, code: str, step: int, window=None) -> np.ndarray:
    """Boolean mask of cells containing `code` at time step `step`."""
    lut = np.array([code in hz for hz in grid.tables[step]], dtype=bool)
    idx = grid.index[step] if window is None else grid.index[step][window]
    idx = np.minimum(idx, len(lut) - 1)
    return lut[idx]


def polygonal(geom):
    """Return only the polygon parts of geom as a MultiPolygon (None if none)."""
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, Polygon):
        return MultiPolygon([geom])
    if isinstance(geom, MultiPolygon):
        return geom
    parts = []
    for g in getattr(geom, "geoms", []):
        sub = polygonal(g)
        if sub is not None:
            parts.extend(sub.geoms)
    return MultiPolygon(parts) if parts else None


def polygonize(mask: np.ndarray, transform) -> list:
    if not mask.any():
        return []
    return [shape(g) for g, v in rasterio.features.shapes(mask.astype(np.uint8), mask=mask, transform=transform) if v == 1]


# --------------------------------------------------------------------------- #
# National WWA file
# --------------------------------------------------------------------------- #
def _to_utc(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, utc=True, errors="coerce", format="ISO8601")


def load_national(shp: Path, by_code, by_name, now: datetime, keep_expired: bool) -> gpd.GeoDataFrame:
    g = gpd.read_file(shp)
    if g.crs is None:
        g = g.set_crs(4326)
    g = g.to_crs(4326)
    g = g[g.geometry.notna() & ~g.geometry.is_empty].copy()
    g["geometry"] = g.geometry.make_valid()
    has_code = g["PHENOM"].notna() & g["SIG"].notna()
    g["HAZ_CODE"] = np.where(has_code, g["PHENOM"].fillna("") + "." + g["SIG"].fillna(""), None)
    onset, issued = _to_utc(g.get("ONSET")), _to_utc(g.get("ISSUANCE"))
    ends, expires = _to_utc(g.get("ENDS")), _to_utc(g.get("EXPIRATION"))
    g["_start"] = onset.fillna(issued)
    g["_end"] = ends.fillna(expires)
    if not keep_expired:
        action = g["VTEC"].fillna("").str.slice(3, 6)
        drop = action.isin(["CAN", "EXP", "UPG"]) | (g["_end"].notna() & (g["_end"] < now))
        if drop.any():
            LOG.info("dropping %d cancelled/expired national features", int(drop.sum()))
        g = g[~drop].copy()
    info = [lookup_hazard(c, p, by_code, by_name) for c, p in zip(g["HAZ_CODE"], g.get("PROD_TYPE"))]
    g["HAZ_NAME"] = [i.name for i in info]
    g["PRIORITY"] = [i.priority for i in info]
    g["COLOR"] = [i.color for i in info]
    g["GEOM_SRC"] = "ZONE"
    g["NDFD_FRAC"] = np.nan
    g["SOURCE"] = "WWA"
    LOG.info("national: %d features", len(g))
    return g


# --------------------------------------------------------------------------- #
# Merge
# --------------------------------------------------------------------------- #
def steps_overlapping(grid: NDFDGrid, start, end, tol: timedelta) -> list[int]:
    if pd.isna(start):
        start = grid.times[0]
    if pd.isna(end):
        end = grid.ends[-1]
    s, e = start - tol, end + tol
    # step i covers [times[i], ends[i]); require positive overlap with [s, e)
    return [i for i, (a, b) in enumerate(zip(grid.times, grid.ends)) if a < e and b > s]


def merge(nat: gpd.GeoDataFrame, grid: NDFDGrid, opts) -> gpd.GeoDataFrame:
    tol = timedelta(minutes=opts.time_tolerance_min)
    cell_area = abs(grid.transform.a * grid.transform.e)
    min_part = opts.min_part_cells * cell_area
    ny, nx = grid.index.shape[1:]

    nat = nat.copy()
    grid_box = box(*rasterio.transform.array_bounds(ny, nx, grid.transform))
    # Only CONUS features (entirely inside CONUS_BBOX) are evaluated against NDFD. Everything
    # else (Alaska, Hawaii, Puerto Rico, Guam, ...) keeps its zone geometry untouched; no
    # reprojection round trip, which would also break polygons crossing the antimeridian.
    w, s_, e, n = opts.conus_bbox
    b = nat.geometry.bounds
    cand = (nat["HAZ_CODE"].isin(grid.codes) & (b.minx >= w) & (b.maxx <= e) & (b.miny >= s_) & (b.maxy <= n))
    natp = nat[cand].to_crs(grid.crs)
    # explained[code] -> list of (steps, window slices, zone cell mask) of national features
    explained: dict[str, list] = {}

    for idx in natp.index:
        code = natp.at[idx, "HAZ_CODE"]
        geom = natp.at[idx, "geometry"]
        if not geom.intersects(grid_box):
            continue
        steps = steps_overlapping(grid, natp.at[idx, "_start"], natp.at[idx, "_end"], tol)
        if not steps:
            continue
        w = rasterio.windows.from_bounds(*geom.bounds, transform=grid.transform)
        c0, r0 = max(int(np.floor(w.col_off)) - 1, 0), max(int(np.floor(w.row_off)) - 1, 0)
        c1 = min(int(np.ceil(w.col_off + w.width)) + 1, nx)
        r1 = min(int(np.ceil(w.row_off + w.height)) + 1, ny)
        if c1 <= c0 or r1 <= r0:
            continue
        win = rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)
        sl = win.toslices()
        wtr = rasterio.windows.transform(win, grid.transform)

        # mark cells this feature explains (for NDFD-only detection)
        zone_cells = rasterio.features.rasterize([geom], out_shape=(int(win.height), int(win.width)),
                                                 transform=wtr, all_touched=True, dtype=np.uint8).astype(bool)
        explained.setdefault(code, []).append((steps, sl, zone_cells))

        # NDFD footprint of this hazard code during the feature's valid window
        footprint = np.zeros(zone_cells.shape, dtype=bool)
        for s in steps:
            footprint |= code_mask(grid, code, s, sl)
        # cells whose centre lies in the zone (fall back to touched cells for tiny zones)
        inside = rasterio.features.rasterize([geom], out_shape=zone_cells.shape, transform=wtr,
                                             all_touched=False, dtype=np.uint8).astype(bool)
        if inside.sum() < 1:
            inside = zone_cells
        hits = int((footprint & inside).sum())
        if hits < opts.min_clip_cells:
            continue  # NDFD does not carry this hazard here -> keep zone geometry
        frac = hits / inside.sum()
        nat.at[idx, "NDFD_FRAC"] = round(float(frac), 4)
        if frac >= opts.full_zone_frac:
            nat.at[idx, "GEOM_SRC"] = "ZONE_FULL"  # NDFD covers whole zone: keep exact zone outline
            continue
        # Keep the exact zone outline and cut out the zone cells NDFD has no hazard in.
        # Holes smaller than min_part are ignored as grid noise.
        holes = [p for p in polygonize(inside & ~footprint, wtr) if p.area >= min_part]
        if holes and opts.simplify_m > 0:
            holes = [h.simplify(opts.simplify_m, preserve_topology=True) for h in holes]
        clipped = geom.difference(unary_union(holes)) if holes else geom
        parts = [p for p in getattr(clipped, "geoms", [clipped]) if p.geom_type == "Polygon" and p.area >= min_part]
        if not parts:
            continue
        new_geom = gpd.GeoSeries([unary_union(parts)], crs=grid.crs).to_crs(4326).iloc[0]
        nat.at[idx, "geometry"] = new_geom
        nat.at[idx, "GEOM_SRC"] = "NDFD_CLIP"

    LOG.info("geometry replaced by NDFD for %d features, full-zone matches %d",
             int((nat.GEOM_SRC == "NDFD_CLIP").sum()), int((nat.GEOM_SRC == "ZONE_FULL").sum()))

    new = ndfd_only_features(grid, explained, opts, min_part)
    if len(new):
        if opts.simplify_m > 0:
            new["geometry"] = new.geometry.simplify(opts.simplify_m, preserve_topology=True)
        new = new.to_crs(4326)
        nat = pd.concat([nat, new], ignore_index=True)
    return gpd.GeoDataFrame(nat, geometry="geometry", crs=4326)


def ndfd_only_features(grid: NDFDGrid, explained, opts, min_part) -> gpd.GeoDataFrame:
    """Polygons for NDFD hazard cells/times that no national feature explains."""
    rows = []
    T = len(grid.times)
    for code in sorted(grid.codes):
        stack = np.zeros(grid.index.shape, dtype=bool)
        for s in range(T):
            stack[s] = code_mask(grid, code, s)
        for steps, sl, zone_cells in explained.get(code, []):
            for s in steps:
                stack[s][sl] &= ~zone_cells
        anyt = stack.any(axis=0)
        if not anyt.any():
            continue
        rr, cc = np.nonzero(anyt)
        r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
        sub = stack[:, r0:r1, c0:c1]
        wtr = grid.transform * rasterio.Affine.translation(c0, r0)
        # group cells by identical on/off time signature
        packed = np.packbits(sub, axis=0)                 # (B, h, w)
        flat = packed.reshape(packed.shape[0], -1).T      # (cells, B)
        flat = np.ascontiguousarray(flat).view(np.dtype((np.void, flat.shape[1]))).ravel()
        uniq, labels = np.unique(flat, return_inverse=True)
        labels = labels.reshape(sub.shape[1:])
        runs: dict[tuple[int, int], list[int]] = {}
        for k, u in enumerate(uniq):
            sig = np.unpackbits(np.frombuffer(u.tobytes(), np.uint8))[:T].astype(bool)
            if not sig.any():
                continue
            d = np.diff(np.concatenate([[0], sig.astype(np.int8), [0]]))
            for a, b in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
                runs.setdefault((int(a), int(b) - 1), []).append(k)
        for (a, b), ks in runs.items():
            m = np.isin(labels, ks)
            polys = [p for p in polygonize(m, wtr) if p.area >= min_part]
            if not polys:
                continue
            rows.append({"HAZ_CODE": code, "PHENOM": code[:2], "SIG": code[3],
                         "_start": pd.Timestamp(grid.times[a]), "_end": pd.Timestamp(grid.ends[b]),
                         "geometry": unary_union(polys)})
    g = gpd.GeoDataFrame(rows, geometry="geometry", crs=grid.crs) if rows else gpd.GeoDataFrame(
        columns=["HAZ_CODE", "_start", "_end", "geometry"], geometry="geometry", crs=grid.crs)
    if len(g):
        info = [lookup_hazard(c, None, opts.by_code, opts.by_name) for c in g["HAZ_CODE"]]
        g["HAZ_NAME"] = [i.name for i in info]
        g["PROD_TYPE"] = g["HAZ_NAME"]
        g["PRIORITY"] = [i.priority for i in info]
        g["COLOR"] = [i.color for i in info]
        g["GEOM_SRC"] = "NDFD_ONLY"
        g["SOURCE"] = "NDFD"
    LOG.info("NDFD-only features: %d", len(g))
    return g


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
OUT_COLS = ["HAZ_CODE", "HAZ_NAME", "PRIORITY", "COLOR", "DRAW_ORDR", "VALID_STRT", "VALID_END",
            "GEOM_SRC", "NDFD_FRAC", "SOURCE",
            "CAP_ID", "MSG_TYPE", "PROD_TYPE", "VTEC", "PHENOM", "SIG", "WFO", "EVENT",
            "ISSUANCE", "EXPIRATION", "ONSET", "ENDS", "URL"]


def finalize(g: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    g = g.copy()
    g["geometry"] = [polygonal(x) for x in g.geometry.make_valid()]
    g = g[g.geometry.notna()].copy()
    fmt = lambda s: s.dt.strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: E731
    g["VALID_STRT"] = fmt(pd.to_datetime(g["_start"], utc=True))
    g["VALID_END"] = fmt(pd.to_datetime(g["_end"], utc=True))
    # lowest priority first, highest priority (1) last -> drawn on top
    g = g.sort_values(["PRIORITY", "_start"], ascending=[False, False], kind="stable").reset_index(drop=True)
    g["DRAW_ORDR"] = np.arange(1, len(g) + 1)
    for c in OUT_COLS:
        if c not in g:
            g[c] = None
    return gpd.GeoDataFrame(g[OUT_COLS + ["geometry"]], geometry="geometry", crs=4326)


def write_outputs(g: gpd.GeoDataFrame, outdir: Path, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=outdir) as tmp:
        tmp = Path(tmp)
        gj = tmp / f"{name}.geojson"
        g.to_file(gj, driver="GeoJSON", engine="pyogrio", COORDINATE_PRECISION=5, RFC7946="YES")
        shp = tmp / f"{name}.shp"
        g.to_file(shp, driver="ESRI Shapefile", engine="pyogrio", ENCODING="UTF-8")
        # move into place (each file replaced atomically)
        for f in tmp.iterdir():
            os.replace(f, outdir / f.name)
    LOG.info("wrote %d features to %s/%s.{geojson,shp}", len(g), outdir, name)


# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=str(HERE / "output"), help="output directory")
    p.add_argument("--workdir", default=str(HERE / "work"), help="download/extract directory")
    p.add_argument("--name", default="merged_hazards", help="output base file name")
    p.add_argument("--priority-csv", default=str(DEFAULT_PRIORITY_CSV))
    p.add_argument("--no-download", action="store_true", help="reuse files already in --workdir")
    p.add_argument("--keep-expired", action="store_true", help="keep CAN/EXP/UPG and expired national features")
    p.add_argument("--time-tolerance-min", type=float, default=60,
                   help="slack (minutes) when comparing national valid times with NDFD steps (default 60)")
    p.add_argument("--min-part-cells", type=float, default=2,
                   help="drop polygon parts smaller than this many NDFD cells (~6.45 km2 each; default 2)")
    p.add_argument("--min-clip-cells", type=float, default=4,
                   help="minimum clipped area (cells) for an NDFD match to replace a zone geometry (default 4)")
    p.add_argument("--full-zone-frac", type=float, default=0.95,
                   help="if NDFD covers at least this fraction of a zone, keep the exact zone outline (default 0.95)")
    p.add_argument("--simplify-m", type=float, default=500,
                   help="geometry simplification tolerance in metres, 0 to disable (default 500)")
    p.add_argument("--conus-bbox", type=float, nargs=4, default=[-135.0, 20.0, -60.0, 53.0],
                   metavar=("W", "S", "E", "N"),
                   help="only features entirely inside this lon/lat box are checked against NDFD "
                        "(default -135 20 -60 53 = CONUS + adjacent waters)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    opts = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if opts.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%SZ")
    logging.Formatter.converter = time.gmtime
    now = datetime.now(timezone.utc)
    work = Path(opts.workdir)
    tgz = work / "current_all.tar.gz"
    grbs = [work / f"ds.wwa.{k}.bin" for k in NDFD_URLS]
    if not opts.no_download:
        download(WWA_URL, tgz)
        for (k, url), dest in zip(NDFD_URLS.items(), grbs):
            download(url, dest)
    shpdir = work / "wwa"
    if shpdir.exists():
        shutil.rmtree(shpdir)
    shp = extract_shapefile(tgz, shpdir)

    opts.by_code, opts.by_name = load_priority_table(Path(opts.priority_csv))
    nat = load_national(shp, opts.by_code, opts.by_name, now, opts.keep_expired)
    grid = read_ndfd(grbs)
    merged = finalize(merge(nat, grid, opts))
    write_outputs(merged, Path(opts.outdir), opts.name)
    unknown = sorted(set(merged.loc[merged.PRIORITY == UNKNOWN_PRIORITY, "HAZ_CODE"].dropna()))
    if unknown:
        LOG.warning("hazards missing from priority table (priority %d): %s", UNKNOWN_PRIORITY, unknown)
    return 0


if __name__ == "__main__":
    sys.exit(main())
