"""Placebo: split each reference (and the pooled reference) at random into two halves of 1 km blocks and treat half A as 'permit', half B as 'reference'.
Same bins, joint valid mask and weighting as compare_big.py. The spread of A-minus-B shows how large an 'excess' arises from block sampling alone."""
import numpy as np, rasterio, geopandas as gpd, json
from rasterio.features import rasterize
from obia_change import idx
from common import *
BS=100; DRAWS=300; rng=np.random.default_rng(5)
P=rasterio.open(ROOT/"data"/"big_pre.tif"); Q=rasterio.open(ROOT/"data"/"big_post.tif"); H,W=P.height,P.width
pre=P.read(); post=Q.read()
ndp,ndq=idx(pre)[0],idx(post)[0]; dN=ndq-ndp; dV=post[8]-pre[8]
valid=np.isfinite(ndp)&np.isfinite(ndq)&np.isfinite(pre[8])&np.isfinite(post[8])
rows,cols=np.indices((H,W)); block=(rows//BS)*1000+(cols//BS)
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
def group(geoms):
    m=rasterize([(g,1) for g in geoms],out_shape=(H,W),transform=P.transform,fill=0).astype(bool)&valid
    return dict(pre=ndp[m],dN=dN[m],dV=dV[m],blk=block[m])
def did(t,k,var,wd=0.1,minpx=200):
    e=[];w=[]
    for lo in np.arange(0.2,0.9+1e-9,wd):
        a=(t["pre"]>=lo)&(t["pre"]<lo+wd); b=(k["pre"]>=lo)&(k["pre"]<lo+wd)
        if a.sum()>minpx and b.sum()>minpx: e.append(np.median(t[var][a])-np.median(k[var][b])); w.append(a.sum())
    return float(np.average(e,weights=w)) if w else np.nan
refs={s:group(list(C[C.site==s].geometry)) for s in SITES}; refs["pooled"]=group(list(C.geometry))
res={}
for name,g in refs.items():
    u=np.unique(g["blk"]); out={"dNDVI":[],"dVH":[]}
    for _ in range(DRAWS):
        half=set(rng.choice(u,len(u)//2,replace=False)); inA=np.isin(g["blk"],list(half))
        A={k:g[k][inA] for k in("pre","dN","dV")}; B={k:g[k][~inA] for k in("pre","dN","dV")}
        out["dNDVI"].append(did(A,B,"dN")); out["dVH"].append(did(A,B,"dV"))
    res[name]={v:dict(mean=float(np.nanmean(x)),sd=float(np.nanstd(x)),p025=float(np.nanpercentile(x,2.5)),p975=float(np.nanpercentile(x,97.5))) for v,x in out.items()}
    print(name,{v:(round(r["mean"],3),round(r["p025"],3),round(r["p975"],3)) for v,r in res[name].items()},flush=True)
json.dump(res,open(ROOT/"out"/"placebo.json","w"),indent=1)
