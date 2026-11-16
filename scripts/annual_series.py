"""Annual dry-season series: matched permit-minus-reference change relative to 2018 (2018-NDVI bins), pooled reference. Pre-trend = 2017/2019 vs 2018 offsets."""
import numpy as np, rasterio, geopandas as gpd, json
from rasterio.features import rasterize
from obia_change import idx
from common import *
REPS=300; rng=np.random.default_rng(5); BS=100
P=rasterio.open(ROOT/"data"/"big_pre.tif"); H,W=P.height,P.width; pre=P.read(); ndp=idx(pre)[0]; vhp=pre[8]
Y={2018:(ndp,vhp),2025:(lambda q:(idx(q)[0],q[8]))(rasterio.open(ROOT/"data"/"big_post.tif").read())}
for y in (2017,2019,2020,2021,2022,2023,2024):
    a=rasterio.open(ROOT/"data"/f"annual_{y}.tif").read(); Y[y]=((a[1]-a[0])/(a[1]+a[0]),a[2])
yrs=sorted(Y); valid=np.isfinite(ndp)&np.isfinite(vhp)
for y in yrs: valid&=np.isfinite(Y[y][0])&np.isfinite(Y[y][1])
rows,cols=np.indices((H,W)); block=(rows//BS)*1000+(cols//BS)
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
def grp(geoms):
    m=rasterize([(g,1) for g in geoms],out_shape=(H,W),transform=P.transform,fill=0).astype(bool)&valid
    d=dict(pre=ndp[m],blk=block[m],n=int(m.sum()))
    for y in yrs: d[("N",y)]=Y[y][0][m]-ndp[m]; d[("V",y)]=Y[y][1][m]-vhp[m]
    return d
K=grp(list(C.geometry)); T={s:grp([site_geom(s).geometry.iloc[0]]) for s in SITES}
def did(t,k,key):
    e=[];w=[]
    for lo in np.arange(0.2,0.9+1e-9,0.1):
        a=(t["pre"]>=lo)&(t["pre"]<lo+.1); b=(k["pre"]>=lo)&(k["pre"]<lo+.1)
        if a.sum()>200 and b.sum()>200: e.append(np.median(t[key][a])-np.median(k[key][b])); w.append(a.sum())
    return float(np.average(e,weights=w)) if w else np.nan
def rs(g):
    u,inv=np.unique(g["blk"],return_inverse=True); ix=[np.where(inv==i)[0] for i in range(len(u))]
    sel=np.concatenate([ix[i] for i in rng.integers(0,len(u),len(u))]); return {k:v[sel] for k,v in g.items() if k!="n"}
res={}; lines=[]
for s in SITES:
    bs=[(rs(T[s]),rs(K)) for _ in range(REPS)]
    for v in "NV":
        BB={}
        for y in yrs:
            est=did(T[s],K,(v,y)); b=[did(a,c,(v,y)) for a,c in bs]; BB[y]=np.array(b); ci=np.nanpercentile(b,[2.5,97.5])
            res.setdefault(s,{}).setdefault(v,{})[y]=[est,float(ci[0]),float(ci[1])]
        d=BB[2025]-BB[2017]; dd=np.nanpercentile(d,[2.5,97.5]); e=res[s][v][2025][0]-res[s][v][2017][0]
        res[s][v]["diff_2025_2017"]=[e,float(dd[0]),float(dd[1])]
        lines.append(f"{s} {v} 2025-2017: {e:+.3f} [{dd[0]:+.3f},{dd[1]:+.3f}]")
    for v,nm in (("N","dNDVI"),("V","dVH")):
        lines.append(f"{s} {nm}: "+"  ".join(f"{y}:{res[s][v][y][0]:+.3f}[{res[s][v][y][1]:+.3f},{res[s][v][y][2]:+.3f}]" for y in yrs))
open(ROOT/"out"/"annual_series.txt","w").write("\n".join(lines)); json.dump(res,open(ROOT/"out"/"annual_series.json","w"),indent=1)
print("\n".join(lines)); print("n valid px (pooled ref)",K["n"],{s:T[s]["n"] for s in SITES})
