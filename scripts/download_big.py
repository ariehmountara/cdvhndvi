"""One shared composite window covering all three sites + all control polygons, built in row strips (same dates/method for treated and control)."""
import sys, numpy as np, rasterio, geopandas as gpd
from rasterio.transform import from_origin
from rasterio.windows import Window
import download_epochs as D
from common import *
C=gpd.read_file(ROOT/"data"/"controls.gpkg")
allg=gpd.GeoSeries([site_geom(s).geometry.iloc[0] for s in SITES]+list(C.geometry),crs=CRS)
b=allg.total_bounds; buf=100
X0,Y0=np.floor((b[0]-buf)/10)*10,np.floor((b[1]-buf)/10)*10; X1,Y1=np.ceil((b[2]+buf)/10)*10,np.ceil((b[3]+buf)/10)*10
W,H=int((X1-X0)/10),int((Y1-Y0)/10); STRIP=400
print("big grid",W,H,flush=True)
import os
if os.environ.get('ANNUAL'):  # extra dry-season (Jun-Oct) years for the annual series: python download_big.py y2019 y2020 ...
    D.EPOCHS.update({f"y{y}":(f"{y}-06-01",f"{y}-10-31") for y in range(2017,2026)})
state={}
def sg(name): return gpd.GeoDataFrame(geometry=[state["box"]],crs=CRS)
def gr(name,buf=0): return state["grid"]
D.site_geom=sg; D.grid=gr
from shapely.geometry import box
for ep in (sys.argv[1:] or D.EPOCHS):
    out=ROOT/"data"/f"big_{ep}.tif"
    with rasterio.open(out,"w",driver="GTiff",height=H,width=W,count=9,dtype="float32",crs=CRS,transform=from_origin(X0,Y1,10,10),compress="deflate",nodata=np.nan,tiled=True,blockxsize=256,blockysize=256,BIGTIFF="YES") as d:
        for r0 in range(0,H,STRIP):
            h=min(STRIP,H-r0); ytop=Y1-r0*10; ybot=ytop-h*10
            state["box"]=box(X0,ybot,X1,ytop); state["grid"]=(X0,ybot,X1,ytop,W,h)
            for attempt in range(3):
                try:
                    c,n=D.s2("big",ep); sar=D.s1("big",ep); break
                except Exception as e: print("retry",r0,repr(e)[:100],flush=True)
            d.write(np.concatenate([c,sar]).astype("float32"),window=Window(0,r0,W,h),indexes=list(range(1,10)))
            print(ep,"strip",r0,"of",H,"done",flush=True)
