"""Covariate rasters on the big_pre/big_post grid (EPSG:32749, 10 m) for multivariable matching:
elevation (Copernicus DEM GLO-30), slope (deg, from the 30 m DEM), distance to the sea (km, from Impact Observatory 2018 water connected to the southern edge),
distance to built area (km, Impact Observatory 2018 class 7), forest-estate function (kws: 1 HL, 2 HP, 3 HPT, 0 other).
Output data/covariates.tif with described bands."""
import os, numpy as np, rasterio, planetary_computer as pc, pystac_client, geopandas as gpd
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.features import rasterize
from scipy import ndimage as ndi
import download_annual as A
from common import *
H, W, X0, Y1 = A.H, A.W, A.X0, A.Y1
TR = from_origin(X0, Y1, 10, 10)
cat = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
LL = list(A.LL)

def warp(href, res, rs, w=None, h=None):
    w = w or int(np.ceil(W * 10 / res)); h = h or int(np.ceil(H * 10 / res))
    with rasterio.open(pc.sign(href)) as src:
        with WarpedVRT(src, crs=CRS, transform=from_origin(X0, Y1, res, res), width=w, height=h, resampling=rs) as v:
            return v.read(1).astype("float32")

# --- DEM at 30 m, slope, then upsample to 10 m
dem_items = list(cat.search(collections=["cop-dem-glo-30"], bbox=LL).items())
w30, h30 = int(np.ceil(W / 3)), int(np.ceil(H / 3))
dem30 = np.full((h30, w30), np.nan, "float32")
for it in dem_items:
    a = warp(it.assets["data"].href, 30, Resampling.bilinear, w30, h30)
    m = np.isfinite(a); dem30[m] = a[m]
gy, gx = np.gradient(np.nan_to_num(dem30, nan=0.0), 30.0)
slope30 = np.degrees(np.arctan(np.hypot(gx, gy)))
up = lambda a: np.repeat(np.repeat(a, 3, axis=0), 3, axis=1)[:H, :W]
elev = up(dem30); slope = up(slope30)

# --- Impact Observatory 2018 land cover: sea distance and built-area distance
lc_items = list(cat.search(collections=["io-lulc-annual-v02"], bbox=LL, datetime="2018-01-01/2018-12-31").items())
lc = np.zeros((H, W), "uint8")
for it in lc_items:
    a = warp(it.assets["data"].href, 10, Resampling.nearest).astype("uint8"); m = a > 0; lc[m] = a[m]
water = lc == 1
lab, n = ndi.label(water)
sea_labels = np.unique(lab[-5:, :][lab[-5:, :] > 0])        # water touching the southern (ocean) edge
sea = np.isin(lab, sea_labels)
coast_km = ndi.distance_transform_edt(~sea) * 10 / 1000.0
built = lc == 7
built_km = ndi.distance_transform_edt(~built) * 10 / 1000.0 if built.any() else np.full((H, W), np.nan, "float32")

# --- forest-estate function
kws = gpd.read_parquet(os.path.expanduser(r"C:\Users\Pongo\.claude\skills\data-ps-jawa\assets\kwshutan_simplified.parquet")).to_crs(CRS)
code = {"HL": 1, "HP": 2, "HPT": 3}
sh = [(g, code[f]) for g, f in zip(kws.geometry, kws.F_KWS) if f in code]
kw = rasterize(sh, out_shape=(H, W), transform=TR, fill=0, dtype="uint8").astype("float32")

bands = [elev, slope, coast_km.astype("float32"), built_km.astype("float32"), kw]
desc = ("elev_m", "slope_deg", "coast_km", "built_km", "kws_code")
with rasterio.open(ROOT / "data" / "covariates.tif", "w", driver="GTiff", height=H, width=W, count=5, dtype="float32", crs=CRS, transform=TR, compress="deflate", nodata=np.nan) as d:
    d.write(np.stack(bands).astype("float32")); d.descriptions = desc
print("covariates done; elev range", float(np.nanmin(elev)), float(np.nanmax(elev)), "sea px", int(sea.sum()), "built px", int(built.sum()), flush=True)
