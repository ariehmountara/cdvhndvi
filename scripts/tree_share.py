"""Share of the 'Trees' class (value 2) in the Impact Observatory annual land cover, per polygon group, for one year."""
import sys, numpy as np, rasterio, geopandas as gpd, planetary_computer as pc
from pystac_client import Client
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.features import rasterize
from common import *
yr = sys.argv[1] if len(sys.argv) > 1 else "2018"
P = rasterio.open(ROOT / "data" / "big_pre.tif"); H, W = P.height, P.width
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
items = list(cat.search(collections=["io-lulc-annual-v02"], bbox=[112.4, -8.6, 113.0, -8.2], datetime=f"{yr}-01-01/{yr}-12-31").items())
out = np.zeros((H, W), np.uint8)
for it in [i for i in items if i.id.endswith(f"-{yr}")]:
    with rasterio.open(pc.sign(it.assets["data"].href)) as src, WarpedVRT(src, crs=P.crs, transform=P.transform, width=W, height=H, resampling=Resampling.nearest, nodata=0) as v:
        a = v.read(1); out = np.where(a > 0, a, out)
C = gpd.read_file(ROOT / "data" / "controls.gpkg")
def share(geoms):
    m = rasterize([(g, 1) for g in geoms], out_shape=(H, W), transform=P.transform, fill=0).astype(bool) & (out > 0)
    return round(float((out[m] == 2).mean() * 100), 1), int(m.sum()) // 100
print(yr, "items used", len([i for i in items if i.id.endswith(f"-{yr}")]))
for s in SITES:
    print("site", s, share([site_geom(s).geometry.iloc[0]]), "| reference", share(list(C[C.site == s].geometry)))
