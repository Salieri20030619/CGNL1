from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,hashlib,numpy as np,pandas as pd
import statsmodels.api as sm
from scipy import stats
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/actomyosin_investigation_20261012_v1');O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
def bh(p,m=None):
 p=np.asarray(p,float);m=len(p) if m is None else m;order=np.argsort(p);q=np.empty(len(p))
 q[order]=np.minimum(1,np.minimum.accumulate((p[order]*m/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def design(v,cols,target=True):
 d=pd.DataFrame({'const':np.ones(len(v))},index=v.index)
 for col in cols:
  if col=='study':
   for k in sorted(v.dataset.unique())[1:]:d['study_'+k]=v.dataset.eq(k).astype(float)
  elif col=='therapy':
   for k in sorted(v.therapy.unique())[1:]:d['therapy_'+k]=v.therapy.eq(k).astype(float)
  else:d[col]=v[col].astype(float)
 if target:d['CGNL1']=v.CGNL1
 return d
def fitted(v,cols,outcome,need_df=5):
 need=[x for x in cols if x not in ['study','therapy']]+['CGNL1',outcome]
 v=v.dropna(subset=need).reset_index(drop=True);X=design(v,cols)
 active=[k for k in X if k=='const' or X[k].std()>1e-10];X=X[active]
 if 'CGNL1' not in X or len(v)-X.shape[1]<need_df or np.linalg.matrix_rank(X)<X.shape[1] or v[outcome].std()<1e-10:return None
 fit=sm.OLS(v[outcome],X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['CGNL1']
 row=dict(n=len(v),residual_df=int(fit.df_resid),slope=float(fit.params.CGNL1),beta=float(fit.params.CGNL1*v.CGNL1.std()/v[outcome].std()),P=float(fit.pvalues.CGNL1),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]))
 return row,v,X
def predict_raw(train,test,cols,outcome,target):
 raw=design(train,cols,target);rt=design(test,cols,target).reindex(columns=raw.columns,fill_value=0)
 keep=[c for c in raw if c=='const' or raw[c].std()>1e-10];raw=raw[keep];rt=rt[keep]
 if target and 'CGNL1' not in raw:return None
 if len(train)-raw.shape[1]<3 or np.linalg.matrix_rank(raw)<raw.shape[1]:return None
 return rt.to_numpy()@np.linalg.lstsq(raw,train[outcome],rcond=None)[0]
def crossval(v,cols,outcome):
 v=v.dropna(subset=['CGNL1',outcome]+[x for x in cols if x not in ['study','therapy']]).reset_index(drop=True);rows=[]
 for i in range(len(v)):
  tr=v.drop(index=i);te=v.iloc[[i]]
  b=predict_raw(tr,te,cols,outcome,False);a=predict_raw(tr,te,cols,outcome,True)
  if b is not None and a is not None:rows.append(dict(patient=v.patient.iloc[i],observed=v[outcome].iloc[i],baseline=float(b[0]),augmented=float(a[0])))
 return pd.DataFrame(rows)
def performance(p):
 if p.empty:return dict(heldout_n=0,MSE_baseline=np.nan,MSE_augmented=np.nan,improvement=np.nan,patients_improved=0)
 b=(p.observed-p.baseline)**2;a=(p.observed-p.augmented)**2
 return dict(heldout_n=len(p),MSE_baseline=float(b.mean()),MSE_augmented=float(a.mean()),improvement=float(1-a.mean()/b.mean()),patients_improved=int((a<b).sum()))
