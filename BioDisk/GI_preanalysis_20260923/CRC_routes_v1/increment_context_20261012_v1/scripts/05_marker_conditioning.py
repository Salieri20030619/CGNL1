from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,hashlib
import numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/increment_context_20261012_v1');O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
path=R/'POST_PRIMARY_MARKER_CONDITIONING.json'
if not path.exists():
 path.write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),prior_results_known=True,trigger='External cycling increment CGNL1=20.8%; ordinary marker VIL1=54.6%, both versus same base. Determine whether CGNL1 adds information after ordinary epithelial markers.',
  outcome='cycling',measurement='pooled_common_genes main definition',comparators=S['ordinary_markers']+['all_six_joint'],
  models='Existing M1 baseline plus each of six original markers individually, plus all six jointly; compare addition of CGNL1. All five scopes retained. No outcome-dependent choice among markers.',
  inference='HC3 t,full rank,residual_df>=5; BH seven marker-conditioning models within each scope including unevaluable hypotheses; no P or predictive claim for unqualified full models.',
  validation='LOPO,training-only fitting; full conditional model must meet inference qualification before CV. Post-primary diagnostic,not independent validation.'),indent=2))
 (R/'POST_PRIMARY_MARKER_FREEZE.json').write_text(json.dumps(dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest()),indent=2))
v=pd.read_csv(O/'patient_profiles_all_variants.tsv',sep='\t')
v=v[v.variant.eq('pooled_common_genes')&v.definition.isin(['discovery_frozen','author_epithelial_QC'])].reset_index(drop=True)
scopes=[('discovery23',v[v.dataset.ne('Ogden2025')]),('all38',v)]+[(k,u) for k,u in v.groupby('dataset')]
models=[];pred=[];elig=[]
def design(train,test,cols):
 mu=train[cols].mean();sd=train[cols].std();keep=sd[sd>1e-10].index.tolist()
 X=pd.DataFrame({'const':np.ones(len(train))},index=train.index);T=pd.DataFrame({'const':np.ones(len(test))},index=test.index)
 for c in keep:X[c]=(train[c]-mu[c])/sd[c];T[c]=(test[c]-mu[c])/sd[c]
 for k in sorted(train.dataset.unique())[1:]:X['study_'+k]=train.dataset.eq(k).astype(float);T['study_'+k]=test.dataset.eq(k).astype(float)
 return X,T,sd
for scope,u in scopes:
 u=u.reset_index(drop=True)
 for marker in S['ordinary_markers']+['all_six_joint']:
  cols=S['base_covariates']+(S['ordinary_markers'] if marker=='all_six_joint' else [marker])
  X,T,sd=design(u,u,cols+['CGNL1'])
  df=len(u)-X.shape[1]
  qualified=df>=5 and np.linalg.matrix_rank(X)==X.shape[1] and 'CGNL1' in X
  elig.append(dict(scope=scope,marker=marker,n=len(u),df=df,status='estimated' if qualified else 'insufficient_df_or_rank'))
  if not qualified:
   models.append(dict(scope=scope,marker=marker,n=len(u),df=df,beta=np.nan,P=np.nan,family='conditioning7_'+scope));continue
  fit=sm.OLS(u.cycling,X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['CGNL1']/sd.CGNL1
  models.append(dict(scope=scope,marker=marker,n=len(u),df=df,beta=fit.params.CGNL1/u.cycling.std(),slope=fit.params.CGNL1/sd.CGNL1,CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=fit.pvalues.CGNL1,family='conditioning7_'+scope))
  for i in range(len(u)):
   tr=u.drop(index=i);te=u.iloc[[i]];X,T,_=design(tr,te,cols);A,Z,_=design(tr,te,cols+['CGNL1'])
   if np.linalg.matrix_rank(A)<A.shape[1] or np.linalg.matrix_rank(X)<X.shape[1]:continue
   b=float((T.to_numpy()@np.linalg.lstsq(X,tr.cycling,rcond=None)[0])[0]);a=float((Z.to_numpy()@np.linalg.lstsq(A,tr.cycling,rcond=None)[0])[0])
   pred.append(dict(scope=scope,marker=marker,patient=te.patient.iloc[0],observed=te.cycling.iloc[0],baseline=b,augmented=a))
t=pd.DataFrame(models);t['FDR']=np.nan
for _,ids in t.groupby('family').groups.items():
 ids=[i for i in ids if np.isfinite(t.loc[i,'P'])];p=t.loc[ids,'P'].to_numpy()
 if len(ids):t.loc[ids,'FDR']=multipletests(np.r_[p,np.ones(7-len(p))],method='fdr_bh')[1][:len(p)]
t.to_csv(O/'E6_marker_conditioned_models.tsv',sep='\t',index=False)
p=pd.DataFrame(pred);p.to_csv(O/'E6_marker_conditioned_predictions.tsv',sep='\t',index=False)
perf=[]
for (scope,marker),g in p.groupby(['scope','marker']):
 b=(g.baseline-g.observed)**2;a=(g.augmented-g.observed)**2
 perf.append(dict(scope=scope,marker=marker,n=len(g),MSE_baseline=b.mean(),MSE_augmented=a.mean(),improvement=1-a.mean()/b.mean(),patients_improved=int((a<b).sum())))
pd.DataFrame(perf).to_csv(O/'E6_marker_conditioned_increment.tsv',sep='\t',index=False)
pd.DataFrame(elig).to_csv(O/'E6_model_eligibility.tsv',sep='\t',index=False)
(R/'E6_COMPLETE.json').write_text(json.dumps(dict(models=len(t),estimated=int(t.P.notna().sum()),predictions=len(p),all_pass=True),indent=2))
print(t[(t.scope=='Ogden2025')|(t.marker=='all_six_joint')].to_string(index=False))
print(pd.DataFrame(perf).query("scope == 'Ogden2025' or marker == 'all_six_joint'").to_string(index=False))
print('MARKER_CONDITIONING_COMPLETE',flush=True)
