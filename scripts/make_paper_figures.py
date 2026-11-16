import json, numpy as np, rasterio, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from rasterio.windows import from_bounds
from shapely.geometry import box
from obia_change import idx
from common import *

OUT = ROOT / "paper"
OUT.mkdir(exist_ok=True)
NAMES = {"harapan_pertiwi": "Harapan Pertiwi", "bhakti_alam_lestari": "Bhakti Alam Lestari",
         "bangkit_sejahtera": "Bangkit Sejahtera + Tunggul Wulung"}
plt.rcParams.update({"font.family": "serif", "font.size": 8})
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
sites = {s: site_geom(s) for s in SITES}

# ---------- Figure 1: study area ----------
X0, X1, Yb, Yt = 657500, 700410, 9065690, 9090210
kec = gpd.read_parquet(KEC_PARQUET).to_crs(CRS)
kec = kec[kec.intersects(box(X0 - 3000, Yb - 3000, X1 + 3000, Yt + 3000))]
fig, ax = plt.subplots(figsize=(7.2, 4.6))
kec.plot(ax=ax, fc="#f3f1ec", ec="#8c8c8c", lw=0.5)
for _, r in kec.iterrows():
    p = r.geometry.representative_point()
    if X0 < p.x < X1 and Yb < p.y < Yt:
        ax.text(p.x, p.y, r.WADMKC, fontsize=5.5, color="#666", ha="center", style="italic")
gpd.GeoSeries([C.geometry.union_all()], crs=CRS).plot(ax=ax, fc="#4c78a8", alpha=0.45, ec="none")
for s, g in sites.items():
    g.plot(ax=ax, fc="#d9480f", alpha=0.8, ec="k", lw=0.6)
lab = {"harapan_pertiwi": (900, 6800), "bhakti_alam_lestari": (4200, 5200), "bangkit_sejahtera": (-3500, 3800)}
for s, g in sites.items():
    p = g.geometry.iloc[0].representative_point()
    dx, dy = lab[s]
    ax.annotate(NAMES[s], (p.x, p.y), xytext=(p.x + dx, p.y + dy), fontsize=6.5, weight="bold", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="k"), bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
ax.set_xlim(X0, X1); ax.set_ylim(Yb, Yt); ax.set_aspect("equal")
ax.set_xticks(np.arange(660000, 700001, 10000)); ax.set_yticks(np.arange(9070000, 9090001, 5000))
ax.set_xticklabels([f"{int(v/1000)}" for v in ax.get_xticks()]); ax.set_yticklabels([f"{int(v/1000)}" for v in ax.get_yticks()])
ax.set_xlabel("Easting (km), UTM zone 49S"); ax.set_ylabel("Northing (km)")
sb_x, sb_y = X0 + 1200, Yb + 1200
ax.plot([sb_x, sb_x + 5000], [sb_y, sb_y], "k", lw=2); ax.text(sb_x + 2500, sb_y + 300, "5 km", ha="center", fontsize=7)
ax.annotate("", xy=(X1 - 1500, Yt - 9400), xytext=(X1 - 1500, Yt - 12400),
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
ax.text(X1 - 1500, Yt - 9100, "N", ha="center", va="bottom", fontsize=9, weight="bold")
ax.text(X0 + 15000, Yb + 900, "Indian Ocean", fontsize=8, style="italic", color="#2b6cb0", ha="center")
ax.legend(handles=[Patch(fc="#d9480f", alpha=0.8, ec="k", label="Permit areas (IPHPS to HKm)"),
                   Patch(fc="#4c78a8", alpha=0.45, label="Reference forest (HP/HL/HPT outside KHDPK and PS)"),
                   Patch(fc="#f3f1ec", ec="#8c8c8c", label="Kecamatan")], loc="upper left", fontsize=6.5, frameon=True)
java = gpd.read_parquet(ADMIN_PARQUET)
java = java[~java.WADMPR.str.contains("Bali", case=False, na=False)]
jl = java.dissolve(by="WADMPR")
ins = ax.inset_axes([0.66, 0.74, 0.33, 0.21])
jl.plot(ax=ins, fc="#d9d6cc", ec="#8c8c8c", lw=0.25)
bb = gpd.GeoSeries([box(X0, Yb, X1, Yt)], crs=CRS).to_crs(4326).total_bounds
ins.add_patch(Rectangle((bb[0], bb[1]), bb[2] - bb[0], bb[3] - bb[1], fc="none", ec="#d9480f", lw=1.2))
ins.set_xlim(104.8, 114.8); ins.set_ylim(-9.0, -5.6); ins.set_aspect("equal")
ins.set_xticks([106, 110, 114]); ins.set_yticks([-8, -6])
ins.set_xticklabels(["106°E", "110°E", "114°E"], fontsize=4.5)
ins.set_yticklabels(["8°S", "6°S"], fontsize=4.5); ins.tick_params(length=1.5, pad=1)
ins.set_title("Java, Indonesia", fontsize=5.5, pad=3)
ins.patch.set_facecolor("#dcebf5")
for sp in ins.spines.values(): sp.set_linewidth(0.6)
plt.tight_layout(); plt.savefig(OUT / "fig1_study_area.png", dpi=200); plt.close()

# ---------- Figure 2: NDVI and VH, 2018 and 2025 ----------
P = rasterio.open(ROOT / "data" / "big_pre.tif"); Q = rasterio.open(ROOT / "data" / "big_post.tif")
fig, axes = plt.subplots(3, 4, figsize=(7.4, 7.8))
cols = ["NDVI 2018", "NDVI 2025", "VH 2018 (dB)", "VH 2025 (dB)"]
for r, s in enumerate(SITES):
    g = sites[s].geometry.iloc[0]; b = g.bounds; m = 250
    win = from_bounds(b[0] - m, b[1] - m, b[2] + m, b[3] + m, P.transform).round_offsets().round_lengths()
    a = P.read(window=win); c = Q.read(window=win)
    ext = [P.window_transform(win) * (0, win.height), P.window_transform(win) * (win.width, 0)]
    extent = [ext[0][0], ext[1][0], ext[0][1], ext[1][1]]
    data = [idx(a)[0], idx(c)[0], a[8], c[8]]
    for k in range(4):
        ax = axes[r, k]
        im = ax.imshow(data[k], extent=extent, cmap="YlGn" if k < 2 else "gray",
                       vmin=0.0 if k < 2 else -20, vmax=0.9 if k < 2 else -8, interpolation="nearest")
        sites[s].boundary.plot(ax=ax, color="#d62728", lw=0.8)
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
        ax.set_xticks([]); ax.set_yticks([])
        if r == 0: ax.set_title(cols[k], fontsize=8)
        if k == 0: ax.set_ylabel(NAMES[s], fontsize=7)
        if r == 2 and k in (1, 3):
            pass
cax1 = fig.add_axes([0.12, 0.06, 0.30, 0.012]); cax2 = fig.add_axes([0.58, 0.06, 0.30, 0.012])
fig.colorbar(plt.cm.ScalarMappable(plt.Normalize(0, 0.9), "YlGn"), cax=cax1, orientation="horizontal").set_label("NDVI", fontsize=7)
fig.colorbar(plt.cm.ScalarMappable(plt.Normalize(-20, -8), "gray"), cax=cax2, orientation="horizontal").set_label("VH backscatter (dB)", fontsize=7)
plt.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.13, wspace=0.04, hspace=0.04)
plt.savefig(OUT / "fig2_composites.png", dpi=200); plt.close()

# ---------- Figure 3: change by baseline NDVI bin ----------
res = json.load(open(ROOT / "out" / "results_tables.json"))["bins"]
fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.6), sharex=True)
for j, s in enumerate(SITES):
    rows = np.array(res[s]); x = rows[:, 0] + 0.05
    for i, (lab_, ci, ri) in enumerate([("Change in NDVI", 3, 4), ("Change in VH (dB)", 5, 6)]):
        ax = axes[i, j]
        ax.plot(x, rows[:, ci], "-o", color="#d9480f", ms=3, lw=1.2, label="Permit area")
        ax.plot(x, rows[:, ri], "-s", color="#4c78a8", ms=3, lw=1.2, label="Reference forest")
        ax.axhline(0, color="#999", lw=0.5)
        ax.set_ylabel(lab_ if j == 0 else "")
        if i == 0: ax.set_title(NAMES[s], fontsize=7.5)
        if i == 1: ax.set_xlabel("2018 NDVI (bin centre)")
axes[0, 0].legend(fontsize=6.5, frameon=False)
plt.tight_layout(); plt.savefig(OUT / "fig3_bins.png", dpi=200); plt.close()
# ---------- Figure 4: annual series ----------
ann = json.load(open(ROOT / "out" / "annual_series.json"))
fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.6), sharex=True)
cl = ["#d9480f", "#4c78a8", "#2f9e44"]
for j, s in enumerate(SITES):
    for i, v in enumerate("NV"):
        ax = axes[i, j]; d = ann[s][v]; yy = sorted(int(k) for k in d if k.isdigit())
        est = [d[str(y)][0] for y in yy]; lo = [d[str(y)][1] for y in yy]; hi = [d[str(y)][2] for y in yy]
        ax.fill_between(yy, lo, hi, color="#d9480f", alpha=0.2, lw=0)
        ax.plot(yy, est, "-o", color="#d9480f", ms=3, lw=1.2); ax.axhline(0, color="#999", lw=0.5)
        ax.axvline(2023, color="#666", lw=0.6, ls=":")
        ax.set_ylabel(("NDVI difference" if v == "N" else "VH difference (dB)") if j == 0 else "")
        if i == 0: ax.set_title(NAMES[s], fontsize=7.5)
        if i == 1: ax.set_xlabel("Dry season (year)")
plt.tight_layout(); plt.savefig(OUT / "fig4_annual.png", dpi=200); plt.close()
print("figures written")
