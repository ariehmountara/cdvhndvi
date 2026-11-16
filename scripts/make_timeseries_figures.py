import numpy as np, rasterio, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rasterio.windows import from_bounds
from obia_change import idx
from common import *
OUT=ROOT/"paper"; plt.rcParams.update({"font.family":"serif","font.size":7})
NAMES={"harapan_pertiwi":"Harapan Pertiwi","bhakti_alam_lestari":"Bhakti Alam Lestari","bangkit_sejahtera":"Bangkit Sejahtera + Tunggul Wulung"}
P=rasterio.open(ROOT/"data"/"big_pre.tif"); Q=rasterio.open(ROOT/"data"/"big_post.tif")
src={2018:(P,lambda a:(idx(a)[0],a[8])),2025:(Q,lambda a:(idx(a)[0],a[8]))}
for y in (2017,2019,2020,2021,2022,2023,2024): src[y]=(rasterio.open(ROOT/"data"/f"annual_{y}.tif"),lambda a:((a[1]-a[0])/(a[1]+a[0]),a[2]))
YRS=sorted(src); sites={s:site_geom(s) for s in SITES}
for k,(lab,cm,vmin,vmax,fn) in enumerate([("NDVI","YlGn",0,0.9,"figS1_timeseries_ndvi.png"),("VH backscatter (dB)","gray",-20,-8,"figS2_timeseries_vh.png")]):
    fig,axes=plt.subplots(3,len(YRS),figsize=(7.4,4.3))
    for r,s in enumerate(SITES):
        b=sites[s].geometry.iloc[0].bounds; m=250
        for c,y in enumerate(YRS):
            ds,f=src[y]; win=from_bounds(b[0]-m,b[1]-m,b[2]+m,b[3]+m,ds.transform).round_offsets().round_lengths()
            d=f(ds.read(window=win))[k]; t=ds.window_transform(win)
            ext=[t*(0,win.height),t*(win.width,0)]; extent=[ext[0][0],ext[1][0],ext[0][1],ext[1][1]]
            ax=axes[r,c]; ax.imshow(d,extent=extent,cmap=cm,vmin=vmin,vmax=vmax,interpolation="nearest")
            sites[s].boundary.plot(ax=ax,color="#d62728",lw=0.5); ax.set_xlim(extent[0],extent[1]); ax.set_ylim(extent[2],extent[3]); ax.set_xticks([]); ax.set_yticks([])
            if r==0: ax.set_title(str(y),fontsize=7)
            if c==0: ax.set_ylabel(NAMES[s].replace(" + ","\n+ "),fontsize=6)
    cax=fig.add_axes([0.35,0.07,0.30,0.018]); fig.colorbar(plt.cm.ScalarMappable(plt.Normalize(vmin,vmax),cm),cax=cax,orientation="horizontal").set_label(lab,fontsize=7)
    plt.subplots_adjust(left=0.08,right=0.99,top=0.94,bottom=0.17,wspace=0.03,hspace=0.04); plt.savefig(OUT/fn,dpi=200); plt.close()
print("ok")
