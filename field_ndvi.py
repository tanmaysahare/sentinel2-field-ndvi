"""Field-level NDVI time series from Sentinel-2 L2A Cloud-Optimized GeoTIFFs on AWS.

Reads only the pixels of one field (a GeoJSON polygon) and a ring of land around it, for every
Sentinel-2 pass in a date range, straight from the public COGs (no AWS account, no download of
whole scenes). Each pass gets a field-level cloud verdict from the scene classification layer.

    python field_ndvi.py example/field.geojson 2026-04-01 2026-06-30 out.csv [44QKJ]

The optional last argument fixes the main MGRS tile (default: the tile with the most passes).

Method (the same one behind https://www.varshakrishi.com/open-data):
  * search the Earth Search STAC API (collection sentinel-2-l2a) for every pass over the field,
    whatever its cloud cover; one scene per date, from the main MGRS tile where it has one
  * read B04 (red), B08 (NIR) and SCL for a window around the field, at 10 m
  * remove the BOA offset (-1000) that processing baseline 04.00+ adds to the reflectance values
  * usable pixel = SCL class 4 (vegetation) or 5 (bare soil); the pass is "clear" for the field
    when at least 90% of the field's pixels are usable
  * field NDVI = mean NDVI of the usable field pixels
  * neighbourhood NDVI = median NDVI of usable vegetated pixels (NDVI > 0) within 500 m of the
    field, excluding a 20 m margin around it
  * a clear reading more than 0.10 below the median of its two clear neighbours on each side is
    flagged as suspected haze (thin cloud the SCL missed)
"""
from __future__ import annotations

import csv
import json
import os
import statistics
import sys
import time

import numpy as np
import rasterio
import requests
from rasterio.features import geometry_mask
from rasterio.warp import transform_geom
from rasterio.windows import from_bounds
from shapely.geometry import mapping, shape

# Public bucket: no credentials. Keep GDAL from listing directories and retry flaky reads.
os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "3")

STAC = "https://earth-search.aws.element84.com/v1/search"
VALID_SCL = (4, 5)          # vegetation, not vegetated (bare soil)
MIN_VALID = 0.90            # share of the field that must be usable
RING_M = 500                # neighbourhood radius, metres
MARGIN_M = 20               # gap between the field and its neighbourhood, metres
HAZE_DROP = 0.10            # NDVI below the local median that marks suspected haze


def search(geom: dict, start: str, end: str, max_cloud: float = 101) -> list[dict]:
    """All sentinel-2-l2a items over the polygon, month by month (the API returns at most 100 per page)."""
    items, page = [], {
        "collections": ["sentinel-2-l2a"],
        "intersects": geom,
        "datetime": f"{start}T00:00:00Z/{end}T23:59:59Z",
        "query": {"eo:cloud_cover": {"lt": max_cloud}},
        "limit": 100,
    }
    url = STAC
    while url:
        for attempt in range(4):
            try:
                r = requests.post(url, json=page, timeout=60)
                r.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 3:
                    raise
                time.sleep(3)
        d = r.json()
        items += d["features"]
        nxt = next((l for l in d.get("links", []) if l.get("rel") == "next"), None)
        url, page = (nxt["href"], nxt.get("body", page)) if nxt else (None, None)
    return items


def one_per_date(items: list[dict], tile: str | None = None) -> list[dict]:
    """One scene per pass date. A field near a tile edge appears in two overlapping MGRS tiles: use
    the main tile (the most common one, or `tile`) and add passes that only the other tile recorded.
    Within a tile, keep the least cloudy scene of the day."""
    tiles = [i["id"].split("_")[1] for i in items]          # e.g. S2B_44QKJ_20261012_0_L2A -> 44QKJ
    main = tile or (max(set(tiles), key=tiles.count) if tiles else None)
    best = {}
    for it, t in zip(items, tiles):
        day, cloud = it["properties"]["datetime"][:10], it["properties"]["eo:cloud_cover"]
        cur = best.get(day)
        rank = (t != main, cloud)                               # main tile first, then less cloud
        if cur is None or rank < cur[0]:
            best[day] = (rank, it)
    return [best[d][1] for d in sorted(best)]


def boa_offset(item: dict) -> float:
    """Baseline 04.00+ L2A values carry a +1000 offset unless Earth Search already removed it."""
    p = item["properties"]
    if p.get("earthsearch:boa_offset_applied") is True:
        return 0.0
    return 1000.0 if str(p.get("s2:processing_baseline", "0")) >= "04.00" else 0.0


def read_window(href: str, bounds, shape_=None):
    with rasterio.open(href) as src:
        win = from_bounds(*bounds, transform=src.transform)
        kw = {"out_shape": shape_, "resampling": rasterio.enums.Resampling.nearest} if shape_ else {}
        data = src.read(1, window=win, boundless=True, fill_value=0, **kw).astype("float32")
        return data, src.window_transform(win) if not shape_ else None


def field_pass(item: dict, field_geom: dict) -> dict:
    A = item["assets"]
    with rasterio.open(A["red"]["href"]) as src:
        crs = src.crs
    field = shape(transform_geom("EPSG:4326", crs, field_geom))          # metres, in the tile's UTM zone
    ring_outer, ring_inner = field.buffer(RING_M), field.buffer(MARGIN_M)
    bounds = ring_outer.bounds

    off = boa_offset(item)
    red, transform = read_window(A["red"]["href"], bounds)
    nir, _ = read_window(A["nir"]["href"], bounds)
    scl, _ = read_window(A["scl"]["href"], bounds, shape_=red.shape)      # SCL is 20 m: resample to 10 m
    red, nir = (red - off) / 10000.0, (nir - off) / 10000.0
    ndvi = (nir - red) / np.maximum(nir + red, 1e-6)
    usable = np.isin(scl, VALID_SCL)

    def inside(g):
        return geometry_mask([mapping(g)], out_shape=red.shape, transform=transform, invert=True)

    in_field = inside(field)
    in_ring = inside(ring_outer) & ~inside(ring_inner)
    valid = float(usable[in_field].mean()) if in_field.any() else 0.0
    rec = {
        "date": item["properties"]["datetime"][:10],
        "scene_id": item["id"],
        "scene_cloud_pct": round(item["properties"]["eo:cloud_cover"], 1),
        "field_pixels": int(in_field.sum()),
        "field_valid_fraction": round(valid, 3),
        "field_ndvi_mean": "",
        "neighbourhood_ndvi_median": "",
    }
    if valid >= MIN_VALID:
        rec["field_ndvi_mean"] = round(float(ndvi[in_field & usable].mean()), 3)
        nb = in_ring & usable & (ndvi > 0)
        if nb.sum() > 500:
            rec["neighbourhood_ndvi_median"] = round(float(np.median(ndvi[nb])), 3)
    return rec


def flag_haze(rows: list[dict]) -> None:
    clear = [r for r in rows if r["field_ndvi_mean"] != ""]
    for i, r in enumerate(clear):
        nb = [clear[j]["field_ndvi_mean"] for j in range(max(0, i - 2), min(len(clear), i + 3)) if j != i]
        r["status"] = "excluded_suspected_haze" if nb and r["field_ndvi_mean"] < statistics.median(nb) - HAZE_DROP else "clear"
    for r in rows:
        r.setdefault("status", "field_not_clear")


def main(geojson: str, start: str, end: str, out: str, tile: str | None = None) -> None:
    g = json.load(open(geojson, encoding="utf-8"))
    field_geom = g["features"][0]["geometry"] if g.get("type") == "FeatureCollection" else g.get("geometry", g)
    items = one_per_date(search(field_geom, start, end), tile)
    print(f"{len(items)} passes over the field between {start} and {end}")
    rows = []
    for it in items:
        try:
            rows.append(field_pass(it, field_geom))
        except rasterio.errors.RasterioIOError as e:
            print("  skipped", it["id"], e)
            continue
        r = rows[-1]
        print(f"  {r['date']}  scene cloud {r['scene_cloud_pct']:5.1f}%  field usable {r['field_valid_fraction']:.0%}  NDVI {r['field_ndvi_mean']}")
    flag_haze(rows)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_clear = sum(r["status"] == "clear" for r in rows)
    print(f"wrote {out}: {len(rows)} passes, {n_clear} clear for the field")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
