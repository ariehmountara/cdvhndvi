"""Imagery-condition figures: (1) study-area cloud/shadow fraction of every candidate Sentinel-2 scene by date, with selections A (tile cloud) and B (study-area cloud) marked;
(2) valid-observation count maps for selection A and B and the NDVI composite A vs B, per year where data/sens_<year>.tif exists. Writes paper/figS3_scene_conditions.png and paper/figS4_valid_counts_<year>.png."""
import os, json, numpy as np, rasterio, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
import planetary_computer as pc
import download_annual as A
import download_sens as S
from common import *
plt.rcParams.update({"font.family": "serif", "font.size": 7.5})
OUT = ROOT / "paper"
# ---------- (1) candidate scene timeline
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), sharey=True)
for ax, yr in zip(axes, (2018, 2025)):
    info = json.load(open(ROOT / "out" / f"sens_scenes_{yr}.json"))
    allc = list(A.D.cat.search(collections=["sentinel-2-l2a"], bbox=list(A.LL), datetime=f"{yr}-06-01/{yr}-10-31", query={"eo:cloud_cover": {"lt": 40}}).items())
    allc = [i for i in allc if i.properties.get("s2:mgrs_tile") == "49LFL"]
    fr = [S.study_cloud_fraction(i) for i in allc]
    inA = [i.id in info["A"] for i in allc]; inB = [i.id in info["B"] for i in allc]
    for it, f, a, b in zip(allc, fr, inA, inB):
        d = it.datetime.timetuple().tm_yday
        ax.scatter(d, f * 100, s=34 if (a or b) else 14, c="#999" if not (a or b) else ("#d9480f" if a and not b else "#4c78a8" if b and not a else "#7b3294"),
                   marker="o" if not a else "s", edgecolors="k" if (a or b) else "none", linewidths=0.5)
    ax.set_title(f"{yr} dry season ({len(allc)} candidate scenes)", fontsize=8); ax.set_xlabel("Day of year"); ax.set_ylim(-3, 105)
    ax.axhspan(-3, 15, color="#e8f0e8", zorder=0)
axes[0].set_ylabel("Cloud/shadow/nodata over study area (%)")
from matplotlib.lines import Line2D
axes[1].legend(handles=[Line2D([], [], marker="s", ls="", mfc="#d9480f", mec="k", label="Selection A only (tile cloud)"), Line2D([], [], marker="o", ls="", mfc="#4c78a8", mec="k", label="Selection B only (study-area cloud)"),
                        Line2D([], [], marker="s", ls="", mfc="#7b3294", mec="k", label="Both"), Line2D([], [], marker="o", ls="", mfc="#999", mec="none", label="Not used")], fontsize=6, loc="upper center", bbox_to_anchor=(-0.05, -0.22), ncol=4, frameon=False)
plt.tight_layout(); plt.savefig(OUT / "figS3_scene_conditions.png", dpi=200, bbox_inches="tight"); plt.close()
# ---------- (2) count maps and composites
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
for yr in (2018, 2025):
    f = ROOT / "data" / f"sens_{yr}.tif"
    if not f.exists() or f.stat().st_size < 1e6: print("skip", yr); continue
    ds = rasterio.open(f); D = {d: i + 1 for i, d in enumerate(ds.descriptions)}
    sub = 4   # downsample for plotting
    rd = lambda k: ds.read(D[k], out_shape=(ds.height // sub, ds.width // sub), resampling=Resampling.average)
    cA, cB, c1 = rd("A_count"), rd("B_count"), rd("S1_count")
    nd = lambda b4, b8: (b8 - b4) / (b8 + b4)
    nA = nd(rd("A_rule_B04"), rd("A_rule_B08")); nB = nd(rd("B_rule_B04"), rd("B_rule_B08"))
    fig, ax = plt.subplots(2, 3, figsize=(7.4, 4.6))
    ext = [ds.bounds.left, ds.bounds.right, ds.bounds.bottom, ds.bounds.top]
    panels = [(cA, "S2 valid scenes, selection A", "viridis", 0, 10), (cB, "S2 valid scenes, selection B", "viridis", 0, 10), (c1, "S1 VH valid scenes", "viridis", 0, max(22, np.nanmax(c1))),
              (nA, "NDVI composite, selection A", "YlGn", 0, 0.9), (nB, "NDVI composite, selection B", "YlGn", 0, 0.9), (nB - nA, "NDVI difference B minus A", "RdBu", -0.15, 0.15)]
    for a, (img, t, cm, lo, hi) in zip(ax.ravel(), panels):
        im = a.imshow(img, extent=ext, cmap=cm, vmin=lo, vmax=hi, interpolation="nearest"); a.set_title(t, fontsize=6.5); a.set_xticks([]); a.set_yticks([])
        gpd.GeoSeries([site_geom(s).geometry.iloc[0] for s in SITES], crs=CRS).boundary.plot(ax=a, color="#d62728", lw=0.4)
        a.set_xlim(ext[0], ext[1]); a.set_ylim(ext[2], ext[3]); fig.colorbar(im, ax=a, fraction=0.04, pad=0.01).ax.tick_params(labelsize=5.5)
    fig.suptitle(f"{yr}: valid observations and composites", fontsize=8); plt.tight_layout(); plt.savefig(OUT / f"figS4_valid_counts_{yr}.png", dpi=200); plt.close()
print("figures done")
