"""Control polygons: forest estate (HP/HL/HPT) NOT in KHDPK and NOT in any PS, within R km of each site (buffered exclusion 300 m)."""
import geopandas as gpd, pandas as pd, os
from shapely.ops import unary_union
from common import *
R=15000
ps=gpd.read_parquet(PS_PARQUET).to_crs(CRS)
kws=gpd.read_parquet(KWS_PARQUET).to_crs(CRS); kh=gpd.read_parquet(KHDPK_PARQUET).to_crs(CRS)
forest=kws[kws.F_KWS.isin(["HP","HL","HPT"])]
excl=unary_union(list(ps.buffer(300).geometry)+list(kh.buffer(0).geometry))
out=[]
for n,i in SITES.items():
    c=ps.iloc[i].geometry; win=c.buffer(R)
    f=gpd.clip(forest,win); f=f.geometry.union_all().difference(excl).difference(c.buffer(300))
    g=gpd.GeoDataFrame({"site":[n]},geometry=[f],crs=CRS).explode(index_parts=False)
    g["ha"]=g.area/1e4; g=g[g.ha>=20]
    print(n,"control pieces",len(g),"total ha",round(g.ha.sum()),"largest",round(g.ha.max()) if len(g) else 0)
    out.append(g)
pd.concat(out).to_file(ROOT/"data"/"controls.gpkg",driver="GPKG")
