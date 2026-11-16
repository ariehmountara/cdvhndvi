"""Leave-2025-out regression of the permit-minus-reference offset on rainfall (two predictors) with 95% prediction intervals; saves out/weather_regression.txt."""
import json, numpy as np, pandas as pd
from scipy import stats
from common import *
W=pd.read_csv(ROOT/"out"/"weather_by_year.csv",index_col=0); ann=json.load(open(ROOT/"out"/"annual_series.json")); L=[]
for s in SITES:
    for v,nm in (("N","NDVI"),("V","VH")):
        y=pd.Series({int(k):x[0] for k,x in ann[s][v].items() if k.isdigit()})
        for x in ("ant30_S2","dry_rain"):
            t=W.index.drop(2025); xx=W.loc[t,x].values; yy=y.loc[t].values; n=len(xx)
            r=stats.linregress(xx,yy); x0=W.loc[2025,x]; e=r.intercept+r.slope*x0
            res=yy-(r.intercept+r.slope*xx); se=np.sqrt((res**2).sum()/(n-2)); pi=stats.t.ppf(.975,n-2)*se*np.sqrt(1+1/n+(x0-xx.mean())**2/((xx-xx.mean())**2).sum())
            L.append(f"{s:22s}{nm:5s}{x:9s} R2={r.rvalue**2:.2f} slope_p={r.pvalue:.2f} expected2025={e:+.3f} +/-{pi:.3f} observed2025={y[2025]:+.3f} residual={y[2025]-e:+.3f} explained={e/y[2025]*100:.0f}% x2025={x0:.0f} trainmax={xx.max():.0f}")
open(ROOT/"out"/"weather_regression.txt","w").write("Fit over 2017-2024 (n=8; the 2018 offset is zero by construction), predict 2025, 95% prediction interval\n"+"\n".join(L)); print("\n".join(L))
