# Code for: Baseline Year and Rainfall Affect Apparent Greening in Community-Managed Social Forestry Areas

Analysis code and derived result tables for a Sentinel-1/Sentinel-2 annual time-series study (dry seasons 2017-2025) with rainfall context, built around a matched comparison of three social forestry permit areas in southern Malang, Indonesia (Harapan Pertiwi, Bhakti Alam Lestari, Bangkit Sejahtera), against nearby non-KHDPK, non-PS forest estate (HP/HL/HPT within 15 km). Dry-season (Jun-Oct) composites for 2018 and 2025, an annual series 2017-2025, and NASA POWER rainfall.

Author: Arieh Mountara, Balai Perhutanan Sosial Yogyakarta, Ministry of Forestry, Republic of Indonesia. Licence: MIT (see `LICENSE`). Cite via `CITATION.cff`.

## Requirements
Python 3.13 and the pinned packages in `requirements.txt` (`pip install -r requirements.txt`). Internet access is needed for the imagery and weather downloads. No account or API key is required (Planetary Computer signing is handled by the `planetary-computer` package). The full 2017-2025 download is large (multi-GB rasters) and takes hours.

## Data availability
- Open imagery, downloaded by the scripts from the Microsoft Planetary Computer STAC API: collections `sentinel-2-l2a`, `sentinel-1-rtc` (ascending orbit 54) and `io-lulc-annual-v02` (Impact Observatory land cover, used for tree share).
- Rainfall and weather: NASA POWER daily data, https://power.larc.nasa.gov (point 112.6 E, -8.35 N, 2016-2025), fetched by `scripts/fetch_power.py`.
- Permit boundary polygons are provided by Balai Perhutanan Sosial Yogyakarta and are **not public**; they are not included. Likewise the derived reference-polygon file `data/controls.gpkg` and all rasters under `data/` are not included. The derived result tables reported in the paper are in `results/`.
- Without the permit boundaries the pipeline cannot be re-run end to end; the code documents exactly what was done, and the files in `results/` hold the numbers reported in the paper.

## Configuration
All paths are set in `scripts/common.py`. Defaults are relative to the repository root (`data/`, `out/`, `paper/` are created automatically). Override with environment variables:

| Variable | Meaning | Default |
|---|---|---|
| `OBIA_ROOT` | project root | parent of `scripts/` |
| `PS_LAYER` | permit-boundary layer, any format geopandas reads (GeoPackage, Shapefile, GeoJSON, GeoParquet); **row index defines the three sites** (`SITES` in `common.py`: rows 299, 298, 303) | `data/ps_definitif.gpkg` |
| `KWS_LAYER` | forest-estate function layer (field `F_KWS`; HP/HL/HPT), only for `make_controls.py` | `data/kwshutan.gpkg` |
| `KHDPK_LAYER` | KHDPK layer, only for `make_controls.py` | `data/khdpk.gpkg` |
| `KEC_LAYER`, `ADMIN_LAYER` | kecamatan and kabupaten/kota boundaries, only for Fig. 1 | `data/kecamatan.gpkg`, `data/admin_kabkota.gpkg` |

The former `*_PARQUET` variable names are still honoured. Adapt `SITES` in `common.py` to the row indices (or edit `site_geom`) of your own layer.
| `BS` | bootstrap block size in pixels (100 = 1 km, 200 = 2 km), `compare_big.py` | 100 |
| `REPS` | bootstrap replicates (`compare_big.py` default 1000, `robustness.py` default 300; the paper used 1000) | see left |

Run from the `scripts/` directory (scripts use `from common import *`). bash: `BS=100 REPS=1000 python compare_big.py`; PowerShell: `$env:BS=100; $env:REPS=1000; python compare_big.py`.

## Run order
1. `python make_controls.py` - reference polygons (forest estate outside KHDPK and PS, within 15 km) -> `data/controls.gpkg`.
2. `python download_big.py` - full 2018 and 2025 composites (the two dates of the baseline contrast) on one 10 m grid -> `data/big_pre.tif`, `data/big_post.tif` (9 bands: 7 Sentinel-2, 2 Sentinel-1). It imports `download_epochs.py`.
3. `python download_annual.py 2017 2019 2020 2021 2022 2023 2024` -> `data/annual_<year>.tif`. If strips fail after retries, `repair_annual.py` recomputes them; its `BAD` dictionary is hard-coded to the strips that failed in the original run and must be edited to match your own log.
4. `BS=100 REPS=1000 python compare_big.py` -> `out/comparison_big.json`; then `BS=200 REPS=1000 python compare_big.py` -> `out/comparison_big_bs200.json` (2 km sensitivity).
5. `python placebo.py` -> `out/placebo.json`, `out/placebo.txt`.
6. `python tree_share.py 2018` and `python tree_share.py 2023` -> `out/tree_share.txt` (the released `results/tree_share.txt` contains both years).
7. `REPS=1000 python robustness.py` -> `out/robustness.json`, `out/robustness.txt`.
8. `python results_tables.py` -> `out/results_tables.json`, `out/results_tables.txt`.
9. `python annual_series.py` -> `out/annual_series.json`, `out/annual_series.txt`.
10. `python fetch_power.py` -> `data/power_daily.json`; `python weather_analysis.py` -> `out/weather_by_year.csv`, `out/weather_corr.txt`, `paper/fig5_weather.png`; `python weather_regression.py` -> `out/weather_regression.txt`.
11. Figures: `python make_paper_figures.py`, `python make_timeseries_figures.py` -> `paper/*.png`.

## Scripts
| Script | Role |
|---|---|
| `common.py` | paths, site definitions, CRS (EPSG:32749), epochs |
| `make_controls.py` | builds reference polygons |
| `download_epochs.py` | composite builder (S2 L2A, S1 RTC), imported by the download scripts |
| `download_big.py`, `download_annual.py`, `repair_annual.py` | download composites |
| `obia_change.py` | only `idx()` (spectral indices) is used by the final analysis; its segmentation/classification code is from an earlier approach not used in the paper |
| `compare_big.py`, `placebo.py`, `robustness.py`, `results_tables.py` | matched comparison, block bootstrap, placebo, sensitivity |
| `tree_share.py` | Impact Observatory tree share |
| `annual_series.py`, `fetch_power.py`, `weather_analysis.py`, `weather_regression.py` | annual series and rainfall |
| `make_paper_figures.py`, `make_timeseries_figures.py` | figures |

## Results (`results/`)
`comparison_big*.json`, `compare_big_1000.txt`, `placebo.json/.txt`, `tree_share.txt`, `robustness.json/.txt`, `results_tables.json/.txt`, `annual_series.json/.txt`, `weather_by_year.csv`, `weather_corr.txt`, `weather_regression.txt`. `comparison_big_reps100.json` is a 100-replicate test run, not used in the paper.

## Caveat
An annual series of nine dry seasons and a two-date baseline contrast, three sites, no field validation: the results are an adjusted association, not a causal effect.
