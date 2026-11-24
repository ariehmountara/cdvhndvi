"""Sensitivity composites for 2018 and 2025 on the big_pre/big_post grid.
Selection A = the 10 least cloudy scenes by tile cloud cover (as in the paper); selection B = the 10 scenes with the least cloud/shadow
over the study mask (permit areas + reference polygons, from SCL at 100 m). For each selection, per-pixel median B04/B08 under three
processing-baseline offset rules (rule = 0.1 for baseline >=04.00 or year >=2022, none = 0, half = 0.05) and per-pixel valid-scene counts;
plus the per-pixel count of valid Sentinel-1 VH scenes (ascending orbit 54).
Output data/sens_<year>.tif with described bands. Usage: python download_sens.py 2018 2025"""
import os, sys, time, json
os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff", GDAL_HTTP_MULTIPLEX="YES")
import numpy as np, rasterio, planetary_computer as pc, geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window
from rasterio.enums import Resampling
from concurrent.futures import ThreadPoolExecutor
import download_annual as A
from common import *

H, W, X0, Y1, STRIP = A.H, A.W, A.X0, A.Y1, 300
TR = from_origin(X0, Y1, 10, 10)
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
GEOMS = [site_geom(s).geometry.iloc[0] for s in SITES] + list(C.geometry)

def read(href, x0, ytop, w, h, rs, tries=4):
    for k in range(tries):
        try:
            with rasterio.open(pc.sign(href.split('?')[0])) as src:
                with WarpedVRT(src, crs=CRS, transform=from_origin(x0, ytop, 10, 10), width=w, height=h, resampling=rs) as v:
                    return v.read(1).astype("float32")
        except Exception as e:
            if k == tries - 1: raise
            time.sleep(2 * (k + 1))

def study_cloud_fraction(it):
    """bad (non-clear) fraction of the study mask from SCL at 100 m"""
    w, h = W // 10, H // 10
    tr = from_origin(X0, Y1, 100, 100)
    mask = rasterize([(g, 1) for g in GEOMS], out_shape=(h, w), transform=tr, fill=0).astype(bool)
    with rasterio.open(pc.sign(it.assets["SCL"].href.split('?')[0])) as src:
        with WarpedVRT(src, crs=CRS, transform=tr, width=w, height=h, resampling=Resampling.nearest) as v:
            scl = v.read(1)
    return float(1 - np.isin(scl[mask], [4, 5, 6, 7]).mean())

def offset_for(it, rule):
    base = it.datetime.year >= 2022 or it.properties.get("s2:processing_baseline", "0") >= "04.00"
    if rule == "none": return 0.0
    if rule == "half": return 0.05 if base else 0.0
    return 0.1 if base else 0.0

def process(year):
    s2, s1 = A.items(year)           # s2 = selection A (tile cloud), s1 = orbit-54 ascending scenes
    allc = list(A.D.cat.search(collections=["sentinel-2-l2a"], bbox=list(A.LL), datetime=f"{year}-06-01/{year}-10-31", query={"eo:cloud_cover": {"lt": 40}}).items())
    allc = [i for i in allc if i.properties.get("s2:mgrs_tile") == "49LFL"]
    with ThreadPoolExecutor(6) as ex: fr = list(ex.map(study_cloud_fraction, allc))
    B = [allc[i] for i in np.argsort(fr)[:10]]
    info = dict(year=year, A=[i.id for i in s2], B=[i.id for i in B], overlap=len(set(i.id for i in s2) & set(i.id for i in B)),
                study_bad_fraction_A=[round(fr[allc.index(next(j for j in allc if j.id == i.id))], 3) for i in s2],
                study_bad_fraction_B=[round(float(f), 3) for f in sorted(fr)[:10]])
    print(year, "overlap A/B:", info["overlap"], "bad A:", info["study_bad_fraction_A"], "bad B:", info["study_bad_fraction_B"], flush=True)
    json.dump(info, open(ROOT / "out" / f"sens_scenes_{year}.json", "w"), indent=1)
    U = {i.id: i for i in list(s2) + B}
    ids_A, ids_B = [i.id for i in s2], [i.id for i in B]
    rules = ["rule"] if year < 2022 else ["rule", "none", "half"]
    desc = []
    for sel in ("A", "B"):
        for r in rules: desc += [f"{sel}_{r}_B04", f"{sel}_{r}_B08"]
    desc += ["A_count", "B_count", "S1_count"]
    out = ROOT / "data" / f"sens_{year}.tif"
    with rasterio.open(out, "w", driver="GTiff", height=H, width=W, count=len(desc), dtype="float32", crs=CRS, transform=TR, compress="deflate", nodata=np.nan, tiled=True, blockxsize=256, blockysize=256, BIGTIFF="YES") as d:
        d.descriptions = tuple(desc)
        for r0 in range(0, H, STRIP):
            h = min(STRIP, H - r0); ytop = Y1 - r0 * 10
            def one(it):
                try:
                    scl = read(it.assets["SCL"].href, X0, ytop, W, h, Resampling.nearest); ok = np.isin(scl, [4, 5, 6, 7])
                    b4 = read(it.assets["B04"].href, X0, ytop, W, h, Resampling.bilinear) / 10000.0
                    b8 = read(it.assets["B08"].href, X0, ytop, W, h, Resampling.bilinear) / 10000.0
                    return it.id, ok, b4, b8
                except Exception as e:
                    print("S2 skip", it.id, "row", r0, repr(e)[:80], flush=True); return None
            def one1(it):
                try:
                    a = read(it.assets["vh"].href, X0, ytop, W, h, Resampling.bilinear); return np.isfinite(a) & (a > 0)
                except Exception as e:
                    print("S1 skip", it.id, repr(e)[:80], flush=True); return None
            with ThreadPoolExecutor(8) as ex:
                got = {}
                for r in ex.map(one, list(U.values())):
                    if r is not None: got[r[0]] = r[1:]
                s1ok = [x for x in ex.map(one1, s1) if x is not None]
            bands = []
            for sel, ids in (("A", ids_A), ("B", ids_B)):
                for rule in rules:
                    b4s, b8s = [], []
                    for i in [k for k in ids if k in got]:
                        ok, b4, b8 = got[i]; off = offset_for(U[i], rule)
                        x4 = np.clip(b4 - off, 0, None); x8 = np.clip(b8 - off, 0, None); x4[~ok] = np.nan; x8[~ok] = np.nan
                        b4s.append(x4); b8s.append(x8)
                    import warnings; warnings.simplefilter("ignore")
                    if not b4s: bands += [np.full((h, W), np.nan, 'float32')] * 2
                    else: bands += [np.nanmedian(np.stack(b4s), axis=0), np.nanmedian(np.stack(b8s), axis=0)]
            cntA = np.sum([got[i][0] for i in ids_A if i in got], axis=0).astype("float32"); cntB = np.sum([got[i][0] for i in ids_B if i in got], axis=0).astype("float32")
            cnt1 = np.sum(s1ok, axis=0).astype("float32") if s1ok else np.zeros((h, W), "float32")
            d.write(np.stack(bands + [cntA, cntB, cnt1]).astype("float32"), window=Window(0, r0, W, h))
            print(year, "strip", r0, "of", H, flush=True)
    print(year, "SENS COMPLETE", flush=True)

if __name__ == "__main__":
    for y in sys.argv[1:]: process(int(y))
    print("SENS ALLDONE", flush=True)
