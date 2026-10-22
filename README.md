# sentinel2-field-ndvi

**A 12-month NDVI time series for one farm field, read straight from the Sentinel-2 L2A Cloud-Optimized GeoTIFFs on AWS, with a field-level cloud verdict for every pass.**

> Listed as a tutorial in the usage examples for the [Sentinel-2 L2A COGs dataset on the Registry of Open Data on AWS](https://registry.opendata.aws/sentinel-2-l2a-cogs/) ([awslabs/open-data-registry#3419](https://github.com/awslabs/open-data-registry/pull/3419)).

Sentinel-2 photographs every field on Earth every few days, and the [Sentinel-2 L2A COGs on the Registry of Open Data on AWS](https://registry.opendata.aws/sentinel-2-l2a-cogs/) make it free to read. For a smallholder field of a hectare or two you need about a hundred 10 m pixels per pass, not whole scenes, and an honest answer to one question: *could the satellite actually see this field today?*

This repository shows how, in one short Python file ([`field_ndvi.py`](field_ndvi.py)) and a notebook with real outputs ([`field_ndvi_tutorial.ipynb`](field_ndvi_tutorial.ipynb)):

1. find **every** Sentinel-2 L2A pass over a field polygon with the [Earth Search](https://earth-search.aws.element84.com/v1) STAC API, whatever its cloud cover
2. read only a small window of B04 (red), B08 (near-infrared) and SCL from each COG over HTTP (no AWS account, no scene downloads)
3. remove the processing-baseline 04.00 BOA offset before computing NDVI
4. decide per pass whether the **field** (not the 110 km scene) was clear: at least 90% of the field's pixels are SCL class 4 (vegetation) or 5 (bare soil)
5. compare the field with the fields within 500 m, and flag thin haze the SCL missed
6. handle fields that sit where two MGRS tiles overlap

It is the method behind the open dataset at **[varshakrishi.com/open-data](https://www.varshakrishi.com/open-data)**: a year of NDVI for one field near Katol, Nagpur district, India (CC BY 4.0).

## Run it

```bash
pip install -r requirements.txt
python field_ndvi.py example/field.geojson 2025-10-07 2026-10-12 out.csv 44QKJ
```

Use your own field by replacing `example/field.geojson` with any GeoJSON polygon (EPSG:4326). The last argument (the main MGRS tile) is optional.

Output: one row per pass with `date, scene_id, scene_cloud_pct, field_pixels, field_valid_fraction, field_ndvi_mean, neighbourhood_ndvi_median, status`, where `status` is `clear`, `field_not_clear` or `excluded_suspected_haze`.

## What the example shows

The example polygon is a 1-hectare square (100 pixels) around the published location of the open-data field. It is **not** that field's boundary, which is withheld to protect the landholder's privacy, so its numbers differ slightly from the published file. From 7 Oct 2025 to 12 Oct 2026 ([`example/field_ndvi_2025-2026.csv`](example/field_ndvi_2025-2026.csv)):

- **91 passes, 58 clear for the field.** 90 come from tile 44QKJ; one pass (28 Aug 2026) was recorded only by the overlapping tile 43QHD.
- **Scene cloud cover is a poor guide to one field.** On 21 Dec 2025 the scene was 47.5% cloud and the field was fully visible; on 5 Jan 2026 the scene was 1.6% cloud and none of the field was usable.
- **A field falling faster than its neighbours.** Between 10 and 20 May 2026 the square's NDVI fell from 0.60 to 0.35 (42%), while the fields within 500 m fell 22%.
- **The monsoon blind spot.** From 17 Aug to 12 Oct 2026 there were 15 passes and only 5 clear views of the field.

## Notes and limits

- NDVI is a greenness index, not a diagnosis: a fall can be water stress, harvest, pruning, pests or a new crop.
- The SCL mask misses thin haze. The haze flag (a clear reading more than 0.10 below the median of its two clear neighbours on each side) catches large dips, not small ones.
- Fields smaller than about half a hectare have too few pure pixels once edges are excluded.
- For the monsoon, add Sentinel-1 radar, which sees through cloud.

## Credits and licence

- Data: contains modified Copernicus Sentinel-2 data 2025-2026, from the [Registry of Open Data on AWS](https://registry.opendata.aws/sentinel-2-l2a-cogs/) (Element 84 Earth Search).
- Code: MIT ([LICENSE](LICENSE)). The published field dataset is CC BY 4.0: cite it as "VarshaKrishi open data (https://www.varshakrishi.com/open-data)".
- Written by [VarshaKrishi](https://www.varshakrishi.com), Nagpur, which builds weekly field-level crop advice for Indian farmers.
