"""Rainfall/weather context for the annual series (NASA POWER daily, 0.5 deg cell) vs scene dates, reference-forest levels and the permit-minus-reference offsets."""
import json, numpy as np, pandas as pd, rasterio, geopandas as gpd
from scipy.stats import spearmanr
from rasterio.features import rasterize
from obia_change import idx
from common import *
import download_annual as A
p=json.load(open(ROOT/"data"/"power_daily.json")); df=pd.DataFrame(p); df.index=pd.to_datetime(df.index,format="%Y%m%d"); df=df.replace(-999,np.nan)
YRS=list(range(2017,2026)); rows=[]
for y in YRS:
    s2,s1=A.items(y); d2=sorted(i.datetime.date() for i in s2); d1=sorted(i.datetime.date() for i in s1)
    ds=df[f"{y}-06-01":f"{y}-10-31"]; wet=df[f"{y-1}-11-01":f"{y}-04-30"]
    ant=lambda dates,n:float(np.mean([df.PRECTOTCORR[pd.Timestamp(d)-pd.Timedelta(days=n):pd.Timestamp(d)-pd.Timedelta(days=1)].sum() for d in dates]))
    rows.append(dict(year=y,dry_rain=ds.PRECTOTCORR.sum(),wet_rain_prev=wet.PRECTOTCORR.sum(),rain_days_dry=int((ds.PRECTOTCORR>=1).sum()),jun_oct_T=ds.T2M.mean(),RH=ds.RH2M.mean(),SW=ds.ALLSKY_SFC_SW_DWN.mean(),
        ant30_S2=ant(d2,30),ant3_S1=ant(d1,3),ant7_S1=ant(d1,7),n_s2=len(d2),n_s1=len(d1),s2_first=str(d2[0]) if d2 else "",s2_last=str(d2[-1]) if d2 else ""))
W=pd.DataFrame(rows).set_index("year")
# reference-forest absolute levels per year + permit offsets
P=rasterio.open(ROOT/"data"/"big_pre.tif"); H,Wd=P.height,P.width; pre=P.read(); post=rasterio.open(ROOT/"data"/"big_post.tif").read()
C=gpd.read_file(ROOT/"data"/"controls.gpkg"); m=rasterize([(g,1) for g in C.geometry],out_shape=(H,Wd),transform=P.transform,fill=0).astype(bool)
lv={2018:(idx(pre)[0],pre[8]),2025:(idx(post)[0],post[8])}
for y in (2017,2019,2020,2021,2022,2023,2024):
    a=rasterio.open(ROOT/"data"/f"annual_{y}.tif").read(); lv[y]=((a[1]-a[0])/(a[1]+a[0]),a[2])
W["ref_ndvi"]=[float(np.nanmedian(lv[y][0][m])) for y in W.index]; W["ref_vh"]=[float(np.nanmedian(lv[y][1][m])) for y in W.index]
ann=json.load(open(ROOT/"out"/"annual_series.json"))
for v,n in (("N","off_ndvi"),("V","off_vh")): W[n]=[np.mean([ann[s][v][str(y)][0] for s in SITES]) for y in W.index]
W.round(3).to_csv(ROOT/"out"/"weather_by_year.csv")
print(W.round(2).to_string())
out=[]
for a in ("dry_rain","wet_rain_prev","rain_days_dry","ant30_S2","ant7_S1","ant3_S1","jun_oct_T","RH","SW"):
    for b in ("ref_ndvi","ref_vh","off_ndvi","off_vh"):
        r,pv=spearmanr(W[a],W[b]); out.append(f"{a:14s} vs {b:9s} rho={r:+.2f} p={pv:.2f}")
print("\n".join(out)); open(ROOT/"out"/"weather_corr.txt","w").write(W.round(2).to_string()+"\n\n"+"\n".join(out))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":7.5})
fig,ax=plt.subplots(1,3,figsize=(7.4,2.5))
ax[0].bar(W.index,W.dry_rain,color=["#d9480f" if y in(2018,2025) else "#9db4cc" for y in W.index]); ax[0].set_ylabel("Jun-Oct rainfall (mm)"); ax[0].set_title("Dry-season rainfall",fontsize=8)
ax[1].plot(W.index,W.ref_ndvi,"-o",ms=3,color="#2f9e44"); ax[1].set_ylabel("Reference median NDVI"); ax[1].set_title("Reference forest NDVI",fontsize=8)
ax[2].scatter(W.dry_rain,W.off_ndvi,s=14,color="#d9480f")
for y in W.index: ax[2].annotate(str(y),(W.dry_rain[y],W.off_ndvi[y]),fontsize=6,xytext=(2,2),textcoords="offset points")
ax[2].set_xlabel("Jun-Oct rainfall (mm)"); ax[2].set_ylabel("Permit-minus-reference NDVI"); ax[2].set_title("Offset vs rainfall",fontsize=8)
plt.tight_layout(); plt.savefig(ROOT/"paper"/"fig5_weather.png",dpi=200)
