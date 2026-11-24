"""Sensitivity of the matched NDVI difference (pooled reference) to (a) scene selection by study-area cloud (B vs A) and (b) the processing-baseline offset
applied to 2025 scenes (rule 0.1 / half 0.05 / none 0), and valid-observation counts per pixel. Needs data/sens_2018.tif and data/sens_2025.tif.
Outputs out/sens_results.json and out/sens_results.txt."""
import json, numpy as np, rasterio, geopandas as gpd
from rasterio.features import rasterize
from common import *
REPS = 200; BS = 100; rng = np.random.default_rng(31)
S18 = rasterio.open(ROOT / "data" / "sens_2018.tif"); S25 = rasterio.open(ROOT / "data" / "sens_2025.tif")
d18 = {d: i + 1 for i, d in enumerate(S18.descriptions)}; d25 = {d: i + 1 for i, d in enumerate(S25.descriptions)}
rd = lambda S, D, k: S.read(D[k])
nd = lambda b4, b8: (b8 - b4) / (b8 + b4)
H, W = S18.height, S18.width
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
def mask(geoms): return rasterize([(g, 1) for g in geoms], out_shape=(H, W), transform=S18.transform, fill=0).astype(bool)
MT = {s: mask([site_geom(s).geometry.iloc[0]]) for s in SITES}; MK = mask(list(C.geometry))
rows, cols = np.indices((H, W)); block = (rows // BS) * 1000 + (cols // BS)
def did(pre_t, d_t, pre_k, d_k, minpx=200):
    e = []; w = []
    for lo in np.arange(0.2, 0.9 + 1e-9, 0.1):
        a = (pre_t >= lo) & (pre_t < lo + .1); b = (pre_k >= lo) & (pre_k < lo + .1)
        if a.sum() > minpx and b.sum() > minpx: e.append(np.median(d_t[a]) - np.median(d_k[b])); w.append(a.sum())
    return float(np.average(e, weights=w)) if w else np.nan
def run(label, sel18, sel25, rule25):
    n18 = nd(rd(S18, d18, f"{sel18}_rule_B04"), rd(S18, d18, f"{sel18}_rule_B08"))
    n25 = nd(rd(S25, d25, f"{sel25}_{rule25}_B04"), rd(S25, d25, f"{sel25}_{rule25}_B08"))
    ok = np.isfinite(n18) & np.isfinite(n25); dN = n25 - n18
    out = {"ref_median_2018": float(np.median(n18[MK & ok])), "ref_median_2025": float(np.median(n25[MK & ok])), "ref_dNDVI_median": float(np.median(dN[MK & ok]))}
    kref = MK & ok; pk = (n18[kref], dN[kref], block[kref])
    for s in SITES:
        t = MT[s] & ok; pt = (n18[t], dN[t], block[t])
        e = did(pt[0], pt[1], pk[0], pk[1])
        def rs(g):
            u, inv = np.unique(g[2], return_inverse=True); ix = [np.where(inv == i)[0] for i in range(len(u))]
            sel = np.concatenate([ix[i] for i in rng.integers(0, len(u), len(u))]); return g[0][sel], g[1][sel]
        b = []
        for _ in range(REPS):
            a1, a2 = rs(pt); c1, c2 = rs(pk); b.append(did(a1, a2, c1, c2))
        ci = np.nanpercentile(b, [2.5, 97.5]); out[s] = dict(est=e, lo=float(ci[0]), hi=float(ci[1]), n_px=int(t.sum()))
    return out
variants = {"A_rule (paper)": ("A", "A", "rule"), "B_rule (study-area cloud selection)": ("B", "B", "rule"), "A_half offset": ("A", "A", "half"), "A_no offset": ("A", "A", "none"),
            "B_half offset": ("B", "B", "half"), "B_no offset": ("B", "B", "none")}
res = {k: run(k, *v) for k, v in variants.items()}
# valid-observation counts per pixel
cnt = {}
for yr, S, D in (("2018", S18, d18), ("2025", S25, d25)):
    for sel in ("A", "B", "S1"):
        c = rd(S, D, f"{sel}_count" if sel != "S1" else "S1_count")
        for g, m in (("permit", np.any([MT[s] for s in SITES], axis=0)), ("reference", MK)):
            v = c[m]; cnt[f"{yr}_{sel}_{g}"] = dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), min=float(v.min()), share_lt3=float((v < 3).mean()), share_lt5=float((v < 5).mean()))
res["valid_observations"] = cnt
json.dump(res, open(ROOT / "out" / "sens_results.json", "w"), indent=1)
L = []
for k, r in res.items():
    if k == "valid_observations": continue
    L.append(f"{k:42s} ref NDVI 2018/2025 {r['ref_median_2018']:.3f}/{r['ref_median_2025']:.3f} (d {r['ref_dNDVI_median']:+.3f}) | " + "  ".join(f"{s[:6]} {r[s]['est']:+.3f} [{r[s]['lo']:+.3f},{r[s]['hi']:+.3f}]" for s in SITES))
L.append("valid observations per pixel:")
for k, v in cnt.items(): L.append(f"  {k:20s} median {v['median']:.0f} p10 {v['p10']:.0f} min {v['min']:.0f} <3: {v['share_lt3']*100:.1f}% <5: {v['share_lt5']*100:.1f}%")
open(ROOT / "out" / "sens_results.txt", "w").write("\n".join(L)); print("\n".join(L))
