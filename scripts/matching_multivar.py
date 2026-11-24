"""Multivariable (coarsened-exact) matching of permit pixels to pooled reference pixels with covariate balance.
Strata: 2018 NDVI bin (0.1) x elevation x slope [x distance to sea] [x distance to built area x forest function].
Models: M0 NDVI only (paper), M1 +elevation+slope, M2 +sea distance, M3 +built-area distance and forest function.
Estimate = permit-pixel-weighted mean over strata (>=MINPX pixels in both groups) of median(permit) - median(reference), for dNDVI and dVH (2025-2018).
Balance = standardised mean difference of each covariate, raw vs after reweighting the reference to the permit stratum distribution.
Spatial block bootstrap (1 km blocks), REPS replicates. Outputs out/multivar.json and out/multivar.txt."""
import os, json, numpy as np, pandas as pd, rasterio, geopandas as gpd
from rasterio.features import rasterize
from obia_change import idx
from common import *
REPS = int(os.environ.get("REPS", 200)); BS = 100; MINPX = int(os.environ.get("MINPX", 200)); rng = np.random.default_rng(23)
P = rasterio.open(ROOT / "data" / "big_pre.tif"); Q = rasterio.open(ROOT / "data" / "big_post.tif"); H, W = P.height, P.width
pre = P.read(); post = Q.read(); ndp, ndq = idx(pre)[0], idx(post)[0]; dN = ndq - ndp; dV = post[8] - pre[8]
cov = rasterio.open(ROOT / "data" / "covariates.tif").read(); elev, slope, coast, built, kws = cov
valid = np.isfinite(ndp) & np.isfinite(ndq) & np.isfinite(pre[8]) & np.isfinite(post[8]) & np.isfinite(elev) & np.isfinite(coast) & np.isfinite(built)
rows, cols = np.indices((H, W)); block = (rows // BS) * 1000 + (cols // BS)
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
cut = lambda a, e: np.digitize(a, e)
B_NDVI = np.minimum(np.floor((ndp - 0.2) / 0.1), 6).astype("int16")           # 0.2-0.9 bins; outside range handled below
inr = (ndp >= 0.2) & (ndp < 0.9)
EL = cut(elev, [25, 50, 100, 200]); SL = cut(slope, [5, 15, 30]); CO = cut(coast, [1, 3, 6]); BU = cut(built, [0.5, 2]); KW = kws.astype(int)
STR = {"M0": [B_NDVI], "M1": [B_NDVI, EL, SL], "M2": [B_NDVI, EL, SL, CO], "M3": [B_NDVI, EL, SL, CO, BU, KW]}
def sid(arrs):
    s = np.zeros((H, W), "int64")
    for a in arrs: s = s * 10 + a.astype("int64")
    return s
SID = {m: sid(a) for m, a in STR.items()}
def group(geoms):
    m = rasterize([(g, 1) for g in geoms], out_shape=(H, W), transform=P.transform, fill=0).astype(bool) & valid & inr
    d = dict(dN=dN[m], dV=dV[m], blk=block[m], n=int(m.sum()))
    for k, v in SID.items(): d[k] = v[m]
    for k, v in dict(ndvi=ndp, elev=elev, slope=slope, coast=coast, built=built, hp=(kws == 2).astype("float32"), hl=(kws == 1).astype("float32")).items(): d["c_" + k] = v[m]
    return d
K = group(list(C.geometry)); T = {s: group([site_geom(s).geometry.iloc[0]]) for s in SITES}
def est(t, k, model, key):
    a = pd.DataFrame({"s": t[model], "y": t[key]}).groupby("s").y.agg(["median", "size"]); b = pd.DataFrame({"s": k[model], "y": k[key]}).groupby("s").y.agg(["median", "size"])
    j = a.join(b, lsuffix="_t", rsuffix="_k", how="inner"); j = j[(j.size_t > MINPX) & (j.size_k > MINPX)]
    if len(j) == 0: return np.nan, 0.0
    return float(np.average(j.median_t - j.median_k, weights=j.size_t)), float(j.size_t.sum() / len(t[model]))
def balance(t, k, model):
    a = pd.Series(t[model]).value_counts(); b = pd.Series(k[model]).value_counts()
    ok = a.index[(a > MINPX) & a.index.isin(b.index[b > MINPX])]
    out = {}
    wk = pd.Series(k[model]).map((a[ok] / b[ok]).to_dict()).fillna(0).values          # reweight reference to permit strata
    mt = np.isin(t[model], ok.values)
    for c in ("ndvi", "elev", "slope", "coast", "built", "hp", "hl"):
        xt, xk = t["c_" + c], k["c_" + c]; sd = np.sqrt((xt.var() + xk.var()) / 2) or 1.0
        raw = (xt.mean() - xk.mean()) / sd
        adj = (xt[mt].mean() - np.average(xk, weights=wk)) / sd if wk.sum() > 0 and mt.any() else np.nan
        out[c] = (float(raw), float(adj))
    return out
def resample(g):
    u, inv = np.unique(g["blk"], return_inverse=True); ix = [np.where(inv == i)[0] for i in range(len(u))]
    sel = np.concatenate([ix[i] for i in rng.integers(0, len(u), len(u))]); return {k: v[sel] for k, v in g.items() if k != "n"}
res = {}; lines = []
for s in SITES:
    res[s] = {}
    bs = [(resample(T[s]), resample(K)) for _ in range(REPS)]
    for m in STR:
        r = {}
        for key, nm in (("dN", "dNDVI"), ("dV", "dVH")):
            e, cov_ = est(T[s], K, m, key); b = [est(a, c, m, key)[0] for a, c in bs]; ci = np.nanpercentile(b, [2.5, 97.5])
            r[nm] = dict(est=e, lo=float(ci[0]), hi=float(ci[1]), coverage=cov_)
        r["balance_smd_raw_adj"] = balance(T[s], K, m); res[s][m] = r
        lines.append(f"{s:22s}{m} cov={r['dNDVI']['coverage']*100:3.0f}%  dNDVI {r['dNDVI']['est']:+.3f} [{r['dNDVI']['lo']:+.3f},{r['dNDVI']['hi']:+.3f}]   dVH {r['dVH']['est']:+.2f} [{r['dVH']['lo']:+.2f},{r['dVH']['hi']:+.2f}]   SMD raw->adj: " + " ".join(f"{c}:{v[0]:+.2f}->{v[1]:+.2f}" for c, v in r["balance_smd_raw_adj"].items()))
        print(lines[-1], flush=True)
json.dump(res, open(ROOT / "out" / "multivar.json", "w"), indent=1); open(ROOT / "out" / "multivar.txt", "w").write("\n".join(lines))
print("MULTIVAR DONE", flush=True)
