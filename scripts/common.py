import geopandas as gpd, numpy as np, pathlib
import os
ROOT=pathlib.Path(os.environ.get("OBIA_ROOT") or pathlib.Path(__file__).resolve().parent.parent)
for _d in ("data","out","paper"): (ROOT/_d).mkdir(exist_ok=True)   # data/ rasters, out/ results, paper/ figures
# Non-public / large inputs; override with environment variables (see README).
# Vector inputs (any format geopandas reads: GeoPackage, Shapefile, GeoJSON or GeoParquet). Override with environment variables (see README).
# PS_LAYER: permit-boundary layer (row order defines the SITES indices below)
def _layer(new,old,default): return os.environ.get(new) or os.environ.get(old) or str(ROOT/"data"/default)
PS_LAYER=_layer("PS_LAYER","PS_PARQUET","ps_definitif.gpkg")
KWS_LAYER=_layer("KWS_LAYER","KWS_PARQUET","kwshutan.gpkg")            # forest-estate function layer, field F_KWS (make_controls.py)
KHDPK_LAYER=_layer("KHDPK_LAYER","KHDPK_PARQUET","khdpk.gpkg")           # KHDPK layer (make_controls.py)
KEC_LAYER=_layer("KEC_LAYER","KEC_PARQUET","kecamatan.gpkg")             # kecamatan boundaries (Fig. 1)
ADMIN_LAYER=_layer("ADMIN_LAYER","ADMIN_PARQUET","admin_kabkota.gpkg")   # kabupaten/kota boundaries (Fig. 1 inset)
def read_vector(path):
    return gpd.read_parquet(path) if str(path).lower().endswith((".parquet",".geoparquet")) else gpd.read_file(path)
SITES={"harapan_pertiwi":299,"bhakti_alam_lestari":298,"bangkit_sejahtera":303}
EPOCHS={"pre":("2018-06-01","2018-10-31"),"post":("2025-06-01","2025-10-31")}
CRS="EPSG:32749"
def site_geom(name):
    g=read_vector(PS_LAYER)
    return g.iloc[[SITES[name]]].to_crs(CRS)
def grid(name,buf=100):
    b=site_geom(name).total_bounds
    x0,y0=np.floor((b[0]-buf)/10)*10,np.floor((b[1]-buf)/10)*10
    x1,y1=np.ceil((b[2]+buf)/10)*10,np.ceil((b[3]+buf)/10)*10
    return x0,y0,x1,y1,int((x1-x0)/10),int((y1-y0)/10)
