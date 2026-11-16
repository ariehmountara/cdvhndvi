"""Treated vs reference comparison from the shared composites (big_pre/big_post): same processing for sites and controls.
Matched on 2018 NDVI bins; one consistent joint valid mask (optical+SAR, both epochs); spatial block bootstrap CIs."""
import numpy as np, rasterio, geopandas as gpd, json, sys, os
from rasterio.features import rasterize
from obia_change import idx
from common import *
BS=int(os.environ.get('BS',100))   # block = BS x BS px (100 = 1 km)
REPS=int(os.environ.get('REPS',1000)); rng=np.random.default_rng(11)
P=rasterio.open(ROOT/"data"/"big_pre.tif"); Q=rasterio.open(ROOT/"data"/"big_post.tif"); H,W=P.height,P.width
pre=P.read(); post=Q.read()
ndp,ndq=idx(pre)[0],idx(post)[0]; dN=ndq-ndp; dV=post[8]-pre[8]
valid=np.isfinite(ndp)&np.isfinite(ndq)&np.isfinite(pre[8])&np.isfinite(post[8])
rows,cols=np.indices((H,W)); block=(rows//BS)*1000+(cols//BS)
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
def group(geoms):
    m=rasterize([(g,1) for g in geoms],out_shape=(H,W),transform=P.transform,fill=0).astype(bool)&valid
    return dict(pre=ndp[m],dN=dN[m],dV=dV[m],blk=block[m],n=int(m.sum()))
T={s:group([site_geom(s).geometry.iloc[0]]) for s in SITES}
K={s:group(list(C[C.site==s].geometry)) for s in SITES}
allc=np.zeros((H,W),bool)
K["pooled"]=group(list(C.geometry))
BINS=lambda wd:np.arange(0.2,0.9+1e-9,wd)
def did(t,k,var,wd=0.1,minpx=200):
    e=[];w=[]
    for lo in BINS(wd):
        a=(t["pre"]>=lo)&(t["pre"]<lo+wd); b=(k["pre"]>=lo)&(k["pre"]<lo+wd)
        if a.sum()>minpx and b.sum()>minpx: e.append(np.median(t[var][a])-np.median(k[var][b])); w.append(a.sum())
    return (float(np.average(e,weights=w)) if w else np.nan), (sum(w)/len(t["pre"]) if w else 0)
def resample(g):
    u,inv=np.unique(g["blk"],return_inverse=True); idxs=[np.where(inv==i)[0] for i in range(len(u))]
    pick=rng.integers(0,len(u),len(u)); sel=np.concatenate([idxs[i] for i in pick])
    return {k:g[k][sel] for k in("pre","dN","dV")}
def boot(t,k,var):
    r=[]
    for _ in range(REPS):
        a=did(resample(t),resample(k),var)[0]; r.append(a)
    return np.nanpercentile(r,[2.5,97.5])
out=[]; print(f"{'site':22s}{'reference':22s}{'dNDVI diff [95% CI]':>26s}{'cov':>5s}{'dVH diff dB [95% CI]':>26s}   n_site_ha  sens(0.05/0.2 bins) dNDVI")
for s in SITES:
    for c in list(SITES)+["pooled"]:
        t,k=T[s],K[c]
        a,ca=did(t,k,"dN"); v,cv=did(t,k,"dV"); ci_a=boot(t,k,"dN"); ci_v=boot(t,k,"dV")
        s05=did(t,k,"dN",0.05)[0]; s20=did(t,k,"dN",0.2)[0]
        tag=c+("*" if c==s else "")
        print(f"{s:22s}{tag:22s}{a:+7.3f} [{ci_a[0]:+.3f},{ci_a[1]:+.3f}]{ca*100:5.0f}{v:+7.2f} [{ci_v[0]:+.2f},{ci_v[1]:+.2f}]   {t['n']//100:6d}   {s05:+.3f}/{s20:+.3f}")
        out.append(dict(site=s,ref=c,own=(c==s),dNDVI=a,dNDVI_ci=list(map(float,ci_a)),cov=ca,dVH=v,dVH_ci=list(map(float,ci_v)),sens_bin005=s05,sens_bin02=s20))
print("* = own-site reference")
print("raw medians (site dNDVI,dVH):",{s:(round(float(np.median(T[s]['dN'])),3),round(float(np.median(T[s]['dV'])),2)) for s in SITES},"| reference:",{c:(round(float(np.median(K[c]['dN'])),3),round(float(np.median(K[c]['dV'])),2)) for c in K})
print("baseline 2018 NDVI median site/ref:",{s:round(float(np.median(T[s]['pre'])),3) for s in SITES},{c:round(float(np.median(K[c]['pre'])),3) for c in K})
json.dump(out,open(ROOT/"out"/("comparison_big.json" if BS==100 else f"comparison_big_bs{BS}.json"),"w"),indent=1)
