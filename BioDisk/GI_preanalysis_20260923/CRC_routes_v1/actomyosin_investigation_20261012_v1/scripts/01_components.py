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

base=R.parent/'epithelial_independence_20261012_v1/results'
v=pd.read_csv(base/'A_model_input.tsv',sep='\t');gs=S['genes']
assert np.allclose(v[gs].mean(1),v.actomyosin)
variants={g:[g] for g in gs};variants.update({'minus_'+g:[x for x in gs if x!=g] for g in gs});variants.update(S['submodules']);variants['full']=gs
rows=[];perf=[];elig=[]
for name,genes in variants.items():v['score_'+name]=v[genes].mean(1)
for scope in ['pooled','GSE245552','GSE285991']:
 vv=v.copy() if scope=='pooled' else v[v.dataset.eq(scope)].copy()
 for model,cols in S['models'].items():
  for name in variants:
   fam='component7' if name in gs else 'dropone7' if name.startswith('minus_') else 'submodule3' if name!='full' else 'full1'
   fit=fitted(vv,cols,'score_'+name)
   if fit is None:elig.append(dict(scope=scope,model=model,outcome=name,status='rank_df_variance_failure'));continue
   row,u,X=fit;rows.append(dict(scope=scope,model=model,outcome=name,family=f'{scope}_{model}_{fam}',**row))
   elig.append(dict(scope=scope,model=model,outcome=name,status='estimated',n=row['n']))
   if scope=='pooled' and model=='M1':perf.append(dict(scope=scope,model=model,outcome=name,**performance(crossval(vv,cols,'score_'+name))))
t=pd.DataFrame(rows);t['FDR']=np.nan
for _,ix in t.groupby('family').groups.items():
 f=t.loc[ix[0],'family'];planned=7 if f.endswith('7') else 3 if f.endswith('3') else 1
 t.loc[ix,'FDR']=bh(t.loc[ix,'P'],planned)
t.to_csv(O/'D1_D2_component_models.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(O/'D1_D2_eligibility.tsv',sep='\t',index=False);pd.DataFrame(perf).to_csv(O/'D2_module_predictive_increment.tsv',sep='\t',index=False)
# Exact additive coefficient decomposition on the same M1 design.
f,u,X=fitted(v,S['models']['M1'],'score_full')
co=[]
for g in gs:
 beta=np.linalg.lstsq(X,u[g],rcond=None)[0][-1];co.append(dict(gene=g,component_slope=beta,contribution_to_full_slope=beta/len(gs),signed_effect_fraction=(beta/len(gs))/f['slope']))
pd.DataFrame(co).to_csv(O/'D1_effect_decomposition.tsv',sep='\t',index=False);assert abs(sum(x['contribution_to_full_slope'] for x in co)-f['slope'])<1e-8
inf=sm.OLS(u.score_full,X).fit().get_influence();influence=pd.DataFrame({'patient':u.patient,'leverage':inf.hat_matrix_diag,'Cooks_distance':inf.cooks_distance[0],'target_DFBETA':inf.dfbetas[:,-1]})
influence.to_csv(O/'D3_influence.tsv',sep='\t',index=False)
lo=[]
for patient in v.patient:
 q=v[~v.patient.eq(patient)].copy();fit=fitted(q,S['models']['M1'],'score_full')
 if fit:lo.append(dict(excluded_patient=patient,**fit[0],**performance(crossval(q,S['models']['M1'],'score_full'))))
lo=pd.DataFrame(lo);lo['FDR_exclusion_sensitivity']=bh(lo.P,len(v));lo.to_csv(O/'D3_leave_patient_out.tsv',sep='\t',index=False)
v.to_csv(O/'D1_D3_model_input.tsv',sep='\t',index=False)
(O/'D1_D3_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':len(t),'patient_exclusions':len(lo),'decomposition_verified':True},indent=2))
print('D1_D3_COMPLETE',len(t),flush=True)
