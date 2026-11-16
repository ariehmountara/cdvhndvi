"""Reviewer-requested checks: effective block counts, inward edge buffers on permit and reference polygons, CI precision.
Same estimator as compare_big.py. Output: out/robustness.json, out/robustness.txt"""
import numpy as np, rasterio, geopandas as gpd, json, os
from rasterio.features import rasterize
from obia_change import idx
from common import *
REPS=int(os.environ.get('REPS',300)); rng=np.random.default_rng(23)
P=rasterio.open(ROOT/"data"/"big_pre.tif"); Q=rasterio.open(ROOT/"data"/"big_post.tif"); H,W=P.height,P.width
pre=P.read(); post=Q.read()
ndp,ndq=idx(pre)[0],idx(post)[0]; dN=ndq-ndp; dV=post[8]-pre[8]
valid=np.isfinite(ndp)&np.isfinite(ndq)&np.isfinite(pre[8])&np.isfinite(post[8])
rows,cols=np.indices((H,W))
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
def group(geoms,BS,d=0):
    gg=[g.buffer(-d) if d else g for g in geoms]; gg=[g for g in gg if not g.is_empty]
    m=rasterize([(g,1) for g in gg],out_shape=(H,W),transform=P.transform,fill=0).astype(bool)&valid
    blk=(rows//BS)*1000+(cols//BS)
    return dict(pre=ndp[m],dN=dN[m],dV=dV[m],blk=blk[m],n=int(m.sum()))
def did(t,k,var,wd=0.1,minpx=200):
    e=[];w=[]
    for lo in np.arange(0.2,0.9+1e-9,wd):
        a=(t["pre"]>=lo)&(t["pre"]<lo+wd); b=(k["pre"]>=lo)&(k["pre"]<lo+wd)
        if a.sum()>minpx and b.sum()>minpx: e.append(np.median(t[var][a])-np.median(k[var][b])); w.append(a.sum())
    return float(np.average(e,weights=w)) if w else np.nan
def resample(g):
    u,inv=np.unique(g["blk"],return_inverse=True); idxs=[np.where(inv==i)[0] for i in range(len(u))]
    sel=np.concatenate([idxs[i] for i in rng.integers(0,len(u),len(u))])
    return {k:g[k][sel] for k in("pre","dN","dV")}
def ci(t,k,var): return np.nanpercentile([did(resample(t),resample(k),var) for _ in range(REPS)],[2.5,97.5])
res={'blocks':{},'buffer':[]}; lines=[]
# 1. effective blocks (1 km and 2 km)
for BS in (100,200):
    for s in SITES:
        res['blocks'][f'{s}_{BS}']=int(len(np.unique(group([site_geom(s).geometry.iloc[0]],BS)['blk'])))
        res['blocks'][f'ref_{s}_{BS}']=int(len(np.unique(group(list(C[C.site==s].geometry),BS)['blk'])))
    res['blocks'][f'ref_pooled_{BS}']=int(len(np.unique(group(list(C.geometry),BS)['blk'])))
lines.append('blocks (valid-pixel blocks):'+json.dumps(res['blocks']))
# 2. inward buffers, 1 km blocks, own and pooled references
for d in (0,100,200):
    T={s:group([site_geom(s).geometry.iloc[0]],100,d) for s in SITES}
    K={s:group(list(C[C.site==s].geometry),100,d) for s in SITES}; K['pooled']=group(list(C.geometry),100,d)
    for s in SITES:
        for c in (s,'pooled'):
            a=did(T[s],K[c],'dN'); v=did(T[s],K[c],'dV'); ca=ci(T[s],K[c],'dN'); cv=ci(T[s],K[c],'dV')
            r=dict(buffer_m=d,site=s,ref=c,dNDVI=a,dNDVI_ci=list(map(float,ca)),dVH=v,dVH_ci=list(map(float,cv)),n_site=T[s]['n'])
            res['buffer'].append(r)
            l=f"buf{d:3d} {s:20s}{c:20s} dNDVI {a:+.3f} [{ca[0]:+.4f},{ca[1]:+.4f}]  dVH {v:+.2f} [{cv[0]:+.2f},{cv[1]:+.2f}]  n={T[s]['n']}"
            print(l,flush=True); lines.append(l)
json.dump(res,open(ROOT/'out'/'robustness.json','w'),indent=1); open(ROOT/'out'/'robustness.txt','w').write('\n'.join(lines))
