"""Online Resource 1: Figs S1-S2 (annual NDVI and VH composites) as one PDF."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mi
from matplotlib.backends.backend_pdf import PdfPages
from common import *
P=ROOT/"paper"
caps=[("figS1_timeseries_ndvi.png","Fig. S1 Dry-season (June-October) NDVI composites, 2017-2025, for the three permit areas (rows). Red lines are the permit boundaries. 2018 and 2025 are full composites; 2017 and 2019-2024 are reduced annual composites (ten least-cloudy Sentinel-2 L2A scenes each)."),
      ("figS2_timeseries_vh.png","Fig. S2 Dry-season Sentinel-1 VH backscatter (dB) composites, 2017-2025 (ascending relative orbit 54, median of all scenes), for the three permit areas. Red lines are the permit boundaries."),
      ("figS3_scene_conditions.png","Fig. S3 Cloud, shadow and no-data fraction over the study area (permit areas and reference polygons) of every candidate Sentinel-2 L2A scene (tile 49LFL, tile cloud cover below 40%) in the 2018 and 2025 dry seasons, from the scene classification layer at 100 m. Squares: scenes used in the paper (selection A, by tile cloud cover); circles: scenes with the smallest study-area fraction (selection B). The shaded band marks 15%."),
      ("figS4_valid_counts_2018.png","Fig. S4a 2018: valid observations per pixel (Sentinel-2 for selections A and B, Sentinel-1 VH) and NDVI composites under both selections. Red lines are the permit boundaries."),
      ("figS4_valid_counts_2025.png","Fig. S4b 2025: as Fig. S4a. Under selection A the southern part of the tile, including the permit areas, has a median of 3 valid Sentinel-2 observations per pixel.")]
with PdfPages(P/"Online_Resource_1.pdf") as pdf:
    for f,c in caps:
        fig=plt.figure(figsize=(8.27,11.69)); ax=fig.add_axes([0.05,0.45,0.9,0.45]); ax.imshow(mi.imread(P/f)); ax.axis("off")
        fig.text(0.07,0.40,"\n".join(__import__("textwrap").wrap(c,95)),fontsize=9,va="top",family="serif")
        fig.text(0.07,0.95,"Online Resource 1 - Baseline year and rainfall affect apparent greening in community-managed social forestry areas",fontsize=8,family="serif")
        pdf.savefig(fig); plt.close(fig)
    import pandas as pd, textwrap
    df=pd.read_csv(ROOT/"out"/"scene_inventory.csv")
    df.columns=["Year","S2 cand.","S2 used","Cloud min","Cloud max","Sen2Cor baselines (used scenes)","S1 in extent","S1 used"]
    fig=plt.figure(figsize=(8.27,11.69)); ax=fig.add_axes([0.05,0.55,0.9,0.3]); ax.axis("off")
    tb=ax.table(cellText=df.values.tolist(),colLabels=list(df.columns),loc="upper center",cellLoc="center"); tb.auto_set_font_size(False); tb.set_fontsize(6.5); tb.scale(1,1.6)
    fig.text(0.07,0.50,"\n".join(textwrap.wrap("Table S1 Scene inventory per dry season (1 June-31 October), Sentinel-2 L2A tile 49LFL and Sentinel-1 RTC, from the Microsoft Planetary Computer STAC query of 10 November 2026. S2 cand. = scenes with cloud cover below 40%; S2 used = the ten least cloudy; Cloud = tile cloud cover of the used scenes (%); S1 used = ascending relative orbit 54 scenes.",95)),fontsize=9,va="top",family="serif")
    pdf.savefig(fig); plt.close(fig)
print("ok")
