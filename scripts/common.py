import geopandas as gpd, numpy as np, pathlib
import os
ROOT=pathlib.Path(os.environ.get("OBIA_ROOT") or pathlib.Path(__file__).resolve().parent.parent)
for _d in ("data","out","paper"): (ROOT/_d).mkdir(exist_ok=True)   # data/ rasters, out/ results, paper/ figures
# Non-public / large inputs; override with environment variables (see README).
# PS_PARQUET: permit-boundary layer (row order defines the SITES indices below)
PS_PARQUET=os.environ.get("PS_PARQUET",str(ROOT/"data"/"ps_definitif_simplified.parquet"))
KWS_PARQUET=os.environ.get("KWS_PARQUET",str(ROOT/"data"/"kwshutan_simplified.parquet"))       # forest-estate function layer (make_controls.py)
KHDPK_PARQUET=os.environ.get("KHDPK_PARQUET",str(ROOT/"data"/"khdpk_simplified.parquet"))      # KHDPK layer (make_controls.py)
KEC_PARQUET=os.environ.get("KEC_PARQUET",str(ROOT/"data"/"kecamatan_simplified.parquet"))      # admin boundaries (Fig. 1)
ADMIN_PARQUET=os.environ.get("ADMIN_PARQUET",str(ROOT/"data"/"admin_kabkota_simplified.parquet"))  # admin boundaries (Fig. 1 inset)
SITES={"harapan_pertiwi":299,"bhakti_alam_lestari":298,"bangkit_sejahtera":303}
EPOCHS={"pre":("2018-06-01","2018-10-31"),"post":("2025-06-01","2025-10-31")}
CRS="EPSG:32749"
def site_geom(name):
    g=gpd.read_parquet(PS_PARQUET)
    return g.iloc[[SITES[name]]].to_crs(CRS)
def grid(name,buf=100):
    b=site_geom(name).total_bounds
    x0,y0=np.floor((b[0]-buf)/10)*10,np.floor((b[1]-buf)/10)*10
    x1,y1=np.ceil((b[2]+buf)/10)*10,np.ceil((b[3]+buf)/10)*10
    return x0,y0,x1,y1,int((x1-x0)/10),int((y1-y0)/10)
