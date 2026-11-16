import numpy as np, rasterio, geopandas as gpd, json
from rasterio.features import rasterize
from obia_change import idx
from common import *
P=rasterio.open(ROOT/"data"/"big_pre.tif"); Q=rasterio.open(ROOT/"data"/"big_post.tif"); H,W=P.height,P.width
pre=P.read(); post=Q.read(); ndp,ndq=idx(pre)[0],idx(post)[0]
valid=np.isfinite(ndp)&np.isfinite(ndq)&np.isfinite(pre[8])&np.isfinite(post[8])
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
def grp(g):
    m=rasterize([(x,1) for x in g],out_shape=(H,W),transform=P.transform,fill=0).astype(bool)
    tot=m.sum(); m&=valid
    return dict(m=m,tot=tot)
G={s:grp([site_geom(s).geometry.iloc[0]]) for s in SITES}
for s in SITES: G["ref_"+s]=grp(list(C[C.site==s].geometry))
G["ref_pooled"]=grp(list(C.geometry))
res={}
print(f"{'group':26s}{'ha':>6s}{'NDVI18 p10/50/90':>20s}{'NDVI25 p10/50/90':>20s}{'VH18':>7s}{'VH25':>7s}{'dNDVI':>7s}{'dVH':>6s}{'<0.5(18)':>9s}{'>=0.7(18)':>10s}{'>=0.7(25)':>10s}")
for k,g in G.items():
    m=g["m"]; a=ndp[m]; b=ndq[m]
    r=dict(ha=int(m.sum()//100),ndvi18=[float(x) for x in np.percentile(a,[10,50,90])],ndvi25=[float(x) for x in np.percentile(b,[10,50,90])],
      vh18=float(np.median(pre[8][m])),vh25=float(np.median(post[8][m])),dN=float(np.median(b-a)),dV=float(np.median(post[8][m]-pre[8][m])),
      lt05_18=float((a<0.5).mean()),ge07_18=float((a>=0.7).mean()),ge07_25=float((b>=0.7).mean()))
    res[k]=r
    print(f"{k:26s}{r['ha']:6d}{r['ndvi18'][0]:7.2f}/{r['ndvi18'][1]:.2f}/{r['ndvi18'][2]:.2f}  {r['ndvi25'][0]:7.2f}/{r['ndvi25'][1]:.2f}/{r['ndvi25'][2]:.2f}{r['vh18']:7.2f}{r['vh25']:7.2f}{r['dN']:7.3f}{r['dV']:6.2f}{100*r['lt05_18']:8.0f}%{100*r['ge07_18']:9.0f}%{100*r['ge07_25']:9.0f}%")
# per-bin table for each site vs its own reference
print("\nper-bin (own reference): bin | n_site_ha n_ref_ha | dNDVI site/ref | dVH site/ref")
bins={}
for s in SITES:
    t=G[s]["m"]; r=G["ref_"+s]["m"]
    print("--",s)
    for lo in np.arange(0.2,0.9,0.1):
        a=t&(ndp>=lo)&(ndp<lo+.1); b=r&(ndp>=lo)&(ndp<lo+.1)
        if a.sum()>200 and b.sum()>200:
            row=[round(lo,1),int(a.sum()//100),int(b.sum()//100),float(np.median((ndq-ndp)[a])),float(np.median((ndq-ndp)[b])),float(np.median((post[8]-pre[8])[a])),float(np.median((post[8]-pre[8])[b]))]
            bins.setdefault(s,[]).append(row)
            print(f"{lo:.1f}-{lo+.1:.1f} | {row[1]:5d} {row[2]:5d} | {row[3]:+.3f}/{row[4]:+.3f} | {row[5]:+.2f}/{row[6]:+.2f}")
json.dump(dict(groups=res,bins=bins),open(ROOT/"out"/"results_tables.json","w"),indent=1)
