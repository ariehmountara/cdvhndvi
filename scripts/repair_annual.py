"""Recompute strips whose 3 retries failed (stale data was written) with per-scene tolerance."""
import sys, rasterio
from rasterio.windows import Window
import download_annual as A
BAD={2020:[0,1200,1600,2000,2400],2021:[2000]}
for yr,rows in BAD.items():
    s2,s1=A.items(yr)
    with rasterio.open(A.ROOT/"data"/f"annual_{yr}.tif","r+") as d:
        for r0 in rows:
            h=min(A.STRIP,A.H-r0); ytop=A.Y1-r0*10
            a=A.strip(s2,s1,A.X0,ytop,A.W,h); d.write(a,window=Window(0,r0,A.W,h)); print(yr,"repaired",r0,flush=True)
print("REPAIRDONE",flush=True)
