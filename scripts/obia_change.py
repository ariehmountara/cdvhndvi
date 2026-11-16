"""Segment (SLIC on stacked pre+post indices), per-segment features, RF per epoch (5 classes), change transition matrix.
Labels: NO field data available -> pseudo-labels from conservative spectral rules on segment means, then RF on full feature set.
Reported accuracy is therefore a self-consistency figure, NOT independent validation."""
import sys, json, numpy as np, rasterio, pandas as pd
from rasterio.features import rasterize, shapes
from skimage.segmentation import slic
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import confusion_matrix, cohen_kappa_score, f1_score, accuracy_score
from common import *
CL=["Waterbody","Built/Bare","Open/Agri","Sparse veg","Dense veg"]
def idx(a):
    B,G,R,RE,N,S1,S2=a[:7]; e=1e-6
    ndvi=(N-R)/(N+R+e); ndwi=(G-N)/(G+N+e); bsi=((S1+R)-(N+B))/((S1+R)+(N+B)+e)
    evi2=2.5*(N-R)/(N+2.4*R+1); osavi=1.16*(N-R)/(N+R+0.16)
    return np.stack([ndvi,ndwi,bsi,evi2,osavi])
def seg_stats(feat,lab,n):
    f=feat.reshape(feat.shape[0],-1); l=lab.ravel(); cnt=np.bincount(l,minlength=n)
    out=[]
    for k in f:
        k=np.nan_to_num(k); m=np.bincount(l,k,n)/np.maximum(cnt,1); s=np.sqrt(np.maximum(np.bincount(l,k*k,n)/np.maximum(cnt,1)-m*m,0)); out+= [m,s]
    return np.stack(out,1)
def rules(ndvi,ndwi,bsi):
    y=np.full(len(ndvi),-1)
    y[(ndwi>0.15)&(ndvi<0.1)]=0
    y[(bsi>0.1)&(ndvi<0.2)&(ndwi<0)]=1
    y[(ndvi>=0.2)&(ndvi<0.4)&(bsi<0.1)&(ndwi<0)]=2   # open / agri / grass
    y[(ndvi>=0.4)&(ndvi<0.6)&(ndwi<0)]=3
    y[(ndvi>=0.7)&(ndwi<0)]=4
    return y
def run(name):
    with rasterio.open(ROOT/"data"/f"{name}_pre.tif") as d: pre=d.read(); tr=d.transform; H,W=d.height,d.width
    with rasterio.open(ROOT/"data"/f"{name}_post.tif") as d: post=d.read()
    geom=site_geom(name).geometry.iloc[0]
    inside=rasterize([(geom,1)],out_shape=(H,W),transform=tr,fill=0).astype(bool)
    ip,iq=idx(pre),idx(post)
    # segmentation input: both epochs' NDVI/BSI/NDWI + NIR/SWIR (scaled) -> segments stable across epochs
    seg_in=np.concatenate([ip[:3],iq[:3],pre[[4,5]],post[[4,5]]]);seg_in=np.nan_to_num(seg_in)
    lo,hi=np.percentile(seg_in,[2,98],axis=(1,2),keepdims=True); seg_in=np.clip((seg_in-lo)/(hi-lo+1e-9),0,1)
    ns=max(int(inside.sum()/60),50)   # ~60 px (0.6 ha) mean segment
    seg=slic(np.moveaxis(seg_in,0,-1),n_segments=ns,compactness=0.08,channel_axis=-1,start_label=0,enforce_connectivity=True,mask=inside)
    n=seg.max()+1
    segi=seg.copy(); segi[seg<0]=0
    feats={}
    for ep,a,ii in(("pre",pre,ip),("post",post,iq)):
        feats[ep]=np.hstack([seg_stats(np.concatenate([a[:9],ii]),segi,n)])
    dfe=seg_stats(np.concatenate([iq-ip,(post[7:9]-pre[7:9])]),segi,n)  # index deltas + VV/VH delta
    area=np.bincount(segi[inside],minlength=n)
    valid=(area>=10)
    pred={}; rep={}
    for ep in("pre","post"):
        X=np.hstack([feats[ep],dfe]); ii=ip if ep=="pre" else iq
        m=seg_stats(ii,segi,n)  # mean/std per index: cols 0,2,4 are means of ndvi,ndwi,bsi
        y=rules(m[:,0],m[:,2],m[:,4]); ok=valid&(y>=0)&np.isfinite(X).all(1)
        rf=RandomForestClassifier(300,class_weight="balanced",n_jobs=-1,random_state=1,min_samples_leaf=2)
        Xl,yl=X[ok],y[ok]
        classes,cn=np.unique(yl,return_counts=True)
        keep=np.isin(yl,classes[cn>=6]); Xl,yl=Xl[keep],yl[keep]
        cv=cross_val_predict(rf,Xl,yl,cv=StratifiedKFold(min(5,np.bincount(yl)[np.bincount(yl)>0].min()),shuffle=True,random_state=1))
        rf.fit(Xl,yl); p=np.full(n,-1); p[valid&np.isfinite(X).all(1)]=rf.predict(X[valid&np.isfinite(X).all(1)]); pred[ep]=p
        labs=sorted(set(yl)); rep[ep]={"n_train":int(len(yl)),"OA":float(accuracy_score(yl,cv)),"kappa":float(cohen_kappa_score(yl,cv)),
            "F1":{CL[l]:float(f) for l,f in zip(labs,f1_score(yl,cv,labels=labs,average=None))},
            "confusion":confusion_matrix(yl,cv,labels=labs).tolist(),"classes":[CL[l] for l in labs],
            "importance_top":None}
    veg=lambda c:(c>=3)
    c0,c1=pred["pre"],pred["post"]
    ch=np.full(n,-1)
    ok=(c0>=0)&(c1>=0)
    ch[ok&veg(c0)&veg(c1)&(c1>=c0)]=2; ch[ok&(~veg(c0))&(~veg(c1))]=3
    ch[ok&((c1>c0)&(veg(c1)))&~(veg(c0)&veg(c1)&(c1==c0))]=0   # gain
    ch[ok&(c1<c0)&veg(c0)]=1                                  # loss
    # SAR corroboration: mean VH delta per segment (col layout of dfe: index deltas 5x(mean,std) then VV, VH)
    dvh=dfe[:,12]; dndvi=dfe[:,0]
    rep["sar_corroborated_gain_ha"]=round(float(np.bincount(segi[inside],weights=None,minlength=n)[(ch==0)&(dvh>0.5)].sum())*0.01,2)
    rep["mean_dNDVI_in_gain"]=round(float(np.nanmean(dndvi[ch==0])),3) if (ch==0).any() else None
    NAMES=["Vegetation gain","Vegetation loss","Stable vegetated","Stable non-vegetated"]
    # pixel-level maps
    mapp=lambda p:np.where(inside,p[segi],-1).astype("int8")
    outs={"class_pre":mapp(c0),"class_post":mapp(c1),"change":mapp(ch),"segments":np.where(inside,seg,-1).astype("int32")}
    for k,v in outs.items():
        with rasterio.open(ROOT/"out"/f"{name}_{k}.tif","w",driver="GTiff",height=H,width=W,count=1,dtype=v.dtype.name,crs=CRS,transform=tr,nodata=-1,compress="deflate") as d: d.write(v,1)
    pa=lambda a:a[inside&(a>=0)]
    tm=pd.crosstab(pd.Series([CL[i] for i in pa(outs["class_pre"])],name="pre"),pd.Series([CL[i] for i in pa(outs["class_post"])],name="post"))*0.01
    tm.round(2).to_csv(ROOT/"out"/f"{name}_transition_ha.csv")
    cs=pd.Series(pa(outs["change"])).map(dict(enumerate(NAMES))).value_counts()*0.01
    rep["change_ha"]=cs.round(2).to_dict(); rep["area_ha"]=float(inside.sum()*0.01); rep["n_segments"]=int(valid.sum())
    rep["class_ha"]={ep:{CL[k]:round(float(v)*.01,2) for k,v in zip(*np.unique(pa(outs["class_"+ep]),return_counts=True))} for ep in("pre","post")}
    json.dump(rep,open(ROOT/"out"/f"{name}_report.json","w"),indent=1)
    print(name,json.dumps({k:rep[k] for k in("area_ha","n_segments","change_ha","class_ha","sar_corroborated_gain_ha","mean_dNDVI_in_gain")}),{e:(round(rep[e]["OA"],3),round(rep[e]["kappa"],3)) for e in("pre","post")})
if __name__=="__main__":
    for s in (sys.argv[1:] or SITES): run(s)
