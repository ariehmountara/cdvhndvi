"""Lean annual dry-season composites (Jun-Oct) on the big_pre/big_post grid: only B04, B08 (NDVI) and Sentinel-1 VH, same scene rules as download_epochs.
Output data/annual_<year>.tif with bands B04, B08, VH_dB. Usage: python download_annual.py 2023 2021 ..."""
import sys, numpy as np, rasterio, geopandas as gpd
from rasterio.transform import from_origin
from rasterio.windows import Window
from rasterio.enums import Resampling
from concurrent.futures import ThreadPoolExecutor
import download_epochs as D
from common import *
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
allg=gpd.GeoSeries([site_geom(s).geometry.iloc[0] for s in SITES]+list(C.geometry),crs=CRS)
b=allg.total_bounds; buf=100
X0,Y0=np.floor((b[0]-buf)/10)*10,np.floor((b[1]-buf)/10)*10; X1,Y1=np.ceil((b[2]+buf)/10)*10,np.ceil((b[3]+buf)/10)*10
W,H=int((X1-X0)/10),int((Y1-Y0)/10); STRIP=400
from shapely.geometry import box
import geopandas as gpd
LL=gpd.GeoSeries([box(X0,Y0,X1,Y1)],crs=CRS).to_crs(4326).total_bounds
def items(year):
    s2=list(D.cat.search(collections=["sentinel-2-l2a"],bbox=list(LL),datetime=f"{year}-06-01/{year}-10-31",query={"eo:cloud_cover":{"lt":40}}).items())
    s2=sorted([i for i in s2 if i.properties.get("s2:mgrs_tile")=="49LFL"],key=lambda i:i.properties["eo:cloud_cover"])[:10]
    s1=[i for i in D.cat.search(collections=["sentinel-1-rtc"],bbox=list(LL),datetime=f"{year}-06-01/{year}-10-31").items()
        if i.properties.get("sat:relative_orbit")==54 and i.properties.get("sat:orbit_state")=="ascending"]
    return s2,s1
def strip(s2,s1,x0,ytop,w,h):
    def one2(it):
      try:
        scl=D.read(it.assets["SCL"].href,x0,ytop,w,h,Resampling.nearest); ok=np.isin(scl,[4,5,6,7])
        if ok.mean()<0.05: return None
        a=np.stack([D.read(it.assets[k].href,x0,ytop,w,h,Resampling.bilinear)/10000.0 for k in("B04","B08")])
        if it.datetime.year>=2022 or it.properties.get("s2:processing_baseline","0")>="04.00": a=a-0.1
        a=np.clip(a,0,None); a[:,~ok]=np.nan; return a
      except Exception as e:
        print('S2 skip',it.id,repr(e)[:80],flush=True); return None
    def one1(it):
        try:
            a=D.read(it.assets["vh"].href,x0,ytop,w,h,Resampling.bilinear); a[a<=0]=np.nan; return a
        except Exception as e: print("S1 err",repr(e)[:100],flush=True); return None
    with ThreadPoolExecutor(8) as ex:
        r2=[r for r in ex.map(one2,s2) if r is not None]; r1=[r for r in ex.map(one1,s1) if r is not None]
    import warnings; warnings.simplefilter("ignore")
    opt=np.nanmedian(np.stack(r2),axis=0) if r2 else np.full((2,h,w),np.nan,"float32")
    vh=10*np.log10(np.nanmedian(np.stack(r1),axis=0)) if r1 else np.full((h,w),np.nan,"float32")
    return np.concatenate([opt,vh[None]]).astype("float32")
if __name__=="__main__":
    print("grid",W,H,flush=True)
    for yr in sys.argv[1:]:
        s2,s1=items(yr); print(yr,"S2",len(s2),"S1",len(s1),flush=True)
        with rasterio.open(ROOT/"data"/f"annual_{yr}.tif","w",driver="GTiff",height=H,width=W,count=3,dtype="float32",crs=CRS,transform=from_origin(X0,Y1,10,10),compress="deflate",nodata=np.nan,tiled=True,blockxsize=256,blockysize=256,BIGTIFF="YES") as d:
            d.descriptions=("B04","B08","VH_dB")
            for r0 in range(0,H,STRIP):
                h=min(STRIP,H-r0); ytop=Y1-r0*10
                for attempt in range(3):
                    try: a=strip(s2,s1,X0,ytop,W,h); break
                    except Exception as e: print("retry",r0,repr(e)[:100],flush=True)
                d.write(a,window=Window(0,r0,W,h)); print(yr,"strip",r0,"of",H,"done",flush=True)
        open(ROOT/"data"/f"annual_{yr}.done","w").write("ok"); print(yr,"COMPLETE",flush=True)
    print("ALLDONE",flush=True)
