"""Download NASA POWER daily data (point 112.6 E, -8.35 N, 2016-2025) -> data/power_daily.json
Format: {parameter: {YYYYMMDD: value}}, as read by weather_analysis.py. No API key needed."""
import json, requests
from common import *
URL="https://power.larc.nasa.gov/api/temporal/daily/point"
par=dict(parameters="PRECTOTCORR,T2M,RH2M,ALLSKY_SFC_SW_DWN",community="AG",longitude=112.6,latitude=-8.35,start="20160101",end="20251231",format="JSON")
r=requests.get(URL,params=par,timeout=300); r.raise_for_status()
d=r.json()["properties"]["parameter"]
json.dump(d,open(ROOT/"data"/"power_daily.json","w"))
print({k:len(v) for k,v in d.items()}, "->", ROOT/"data"/"power_daily.json")
