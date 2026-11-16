"""Scene inventory per dry season from the Planetary Computer STAC: candidates, used scenes, cloud range, processing baselines, S1 counts; saves out/scene_inventory.csv."""
import pandas as pd, collections
import download_annual as A
from common import *
rows=[]
for y in range(2017,2026):
    s2,s1=A.items(y)
    allc=list(A.D.cat.search(collections=["sentinel-2-l2a"],bbox=list(A.LL),datetime=f"{y}-06-01/{y}-10-31",query={"eo:cloud_cover":{"lt":40}}).items())
    allc=[i for i in allc if i.properties.get("s2:mgrs_tile")=="49LFL"]
    s1all=[i for i in A.D.cat.search(collections=["sentinel-1-rtc"],bbox=list(A.LL),datetime=f"{y}-06-01/{y}-10-31").items()]
    cl=[i.properties["eo:cloud_cover"] for i in s2]; bl=collections.Counter(i.properties.get("s2:processing_baseline","?") for i in s2)
    rows.append(dict(year=y,s2_candidates_cloud_lt40=len(allc),s2_used=len(s2),cloud_min=round(min(cl),1),cloud_max=round(max(cl),1),baselines=dict(bl),s1_candidates=len(s1all),s1_used_orbit54_asc=len(s1)))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(ROOT/"out"/"scene_inventory.csv",index=False)
