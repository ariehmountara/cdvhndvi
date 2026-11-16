"""Sentinel-2 L2A + Sentinel-1 RTC dry-season median composites per site/epoch (Planetary Computer, no auth)."""
import os
os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",GDAL_HTTP_MULTIPLEX="YES")
import sys, numpy as np, rasterio, planetary_computer as pc, pystac_client
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from concurrent.futures import ThreadPoolExecutor
from common import *
cat=pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=pc.sign_inplace)
BANDS=["B02","B03","B04","B05","B08","B11","B12"]
def read(href,x0,y1,w,h,rs):
    with rasterio.open(pc.sign(href)) as src:
        with WarpedVRT(src,crs=CRS,transform=from_origin(x0,y1,10,10),width=w,height=h,resampling=rs) as v:
            return v.read(1).astype("float32")
def s2(name,ep):
    x0,y0,x1,y1,w,h=grid(name); ll=site_geom(name).to_crs(4326).total_bounds
    its=list(cat.search(collections=["sentinel-2-l2a"],bbox=list(ll),datetime="/".join(EPOCHS[ep]),
        query={"eo:cloud_cover":{"lt":40}}).items())
    its=[i for i in its if i.properties.get("s2:mgrs_tile")=="49LFL"]  # all sites lie in this tile; avoid neighbour-tile scenes
    its=sorted(its,key=lambda i:i.properties["eo:cloud_cover"])[:10]
    print(name,ep,"S2 scenes",len(its),flush=True)
    def one(it):
        scl=read(it.assets["SCL"].href,x0,y1,w,h,Resampling.nearest)
        ok=np.isin(scl,[4,5,6,7])
        if ok.mean()<0.05: return None
        arr=np.stack([read(it.assets[b].href,x0,y1,w,h,Resampling.bilinear)/10000.0 for b in BANDS])
        if it.datetime.year>=2022 or it.properties.get("s2:processing_baseline","0")>="04.00": arr=arr-0.1  # BOA offset
        arr=np.clip(arr,0,None)
        arr[:,~ok]=np.nan; return arr
    with ThreadPoolExecutor(6) as ex: res=[r for r in ex.map(one,its) if r is not None]
    comp=np.nanmedian(np.stack(res),axis=0); cnt=np.sum([~np.isnan(r[0]) for r in res],axis=0)
    return comp,cnt
def s1(name,ep):
    x0,y0,x1,y1,w,h=grid(name); ll=site_geom(name).to_crs(4326).total_bounds
    its=list(cat.search(collections=["sentinel-1-rtc"],bbox=list(ll),datetime="/".join(EPOCHS[ep])).items())
    its=[i for i in its if i.properties.get("sat:relative_orbit")==54 and i.properties.get("sat:orbit_state")=="ascending"]
    print(name,ep,"S1 scenes",len(its),flush=True)
    if not its: return np.full((2,h,w),np.nan,"float32")
    def one(it):
        try:
            a=np.stack([read(it.assets[p].href,x0,y1,w,h,Resampling.bilinear) for p in("vv","vh")]); a[a<=0]=np.nan; return a
        except Exception as e: print("S1 read err",repr(e)[:120],flush=True); return None
    with ThreadPoolExecutor(6) as ex: res=[r for r in ex.map(one,its) if r is not None]
    a=np.nanmedian(np.stack(res),axis=0)
    return 10*np.log10(a)
if __name__=="__main__":
    for name in (sys.argv[1:] or SITES):
        x0,y0,x1,y1,w,h=grid(name)
        for ep in EPOCHS:
            c,n=s2(name,ep); sar=s1(name,ep)
            with rasterio.open(ROOT/"data"/f"{name}_{ep}.tif","w",driver="GTiff",height=h,width=w,count=9,dtype="float32",
                crs=CRS,transform=from_origin(x0,y1,10,10),compress="deflate",nodata=np.nan) as d:
                d.write(np.concatenate([c,sar]).astype("float32")); d.descriptions=tuple(BANDS+["VV_dB","VH_dB"])
            print(name,ep,"done",c.shape,"min valid obs",int(n.min()),flush=True)
