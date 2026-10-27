from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,hashlib
import numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'increment_context_20261012_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
amend=R/'POST_PRIMARY_STATE_INCREMENT.json'
if not amend.exists():
 a=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reason='Four planned state association tests showed an external cycling nominal signal: P=.016430,FDR=.065721; discovery cycling coefficients negative. Test increment for all four states without selection.',
 prior_results_known=True,independent_replication=False,
 outcomes=S['modifiers'],measurement='pooled_common_genes, author epithelial external main definition',
 baseline='depth, epithelial identity, mature_colonocyte; drop mature_colonocyte only when it is the outcome',
 scope='original two individual cohorts,discovery23,external15,all38; patient leave-one-out and three-cohort leave-one-cohort-out',
 comparators='same frozen six ordinary markers and20 expression-matched genes, no new matching',
 target_FDR='four state outcomes within each scope; benchmark27 predictors within scope/outcome',
 inference='HC3 t,full rank,residual_df>=5. Patient-held-out and cohort-held-out errors descriptive, no nominal prediction P.',
 classification='Post-primary exploratory analysis, same patient data; no independent confirmation or causal claim')
 amend.write_text(json.dumps(a,indent=2))
 (R/'POST_PRIMARY_STATE_FREEZE.json').write_text(json.dumps(dict(sha256=hashlib.sha256(amend.read_bytes()).hexdigest()),indent=2))
assert hashlib.sha256(amend.read_bytes()).hexdigest()==json.loads((R/'POST_PRIMARY_STATE_FREEZE.json').read_text())['sha256']
v=pd.read_csv(O/'patient_profiles_all_variants.tsv',sep='\t')
v=v[v.variant.eq('pooled_common_genes')&v.definition.isin(['discovery_frozen','author_epithelial_QC'])].reset_index(drop=True)
control=json.loads((B/'epithelial_independence_20261012_v1/results/A_CONTROL_FREEZE.json').read_text())['controls']
predictors=['CGNL1']+S['ordinary_markers']+control
def matrices(train,test,cols,predictor=None,study=True):
 cols=list(cols)+([predictor] if predictor else [])
 mu=train[cols].mean();sd=train[cols].std()
 keep=sd[sd>1e-10].index.tolist()
 if predictor and predictor not in keep:return None
 x=(train[keep]-mu[keep])/sd[keep];t=(test[keep]-mu[keep])/sd[keep]
 X=pd.DataFrame({'const':np.ones(len(train))},index=train.index);T=pd.DataFrame({'const':np.ones(len(test))},index=test.index)
 for c in keep:X[c]=x[c];T[c]=t[c]
 if study:
  lev=sorted(train.dataset.unique())
  if not set(test.dataset)<=set(lev):return None
  for k in lev[1:]:X['study_'+k]=train.dataset.eq(k).astype(float);T['study_'+k]=test.dataset.eq(k).astype(float)
 if not np.isfinite(X).all().all() or not np.isfinite(T).all().all() or np.linalg.matrix_rank(X)<X.shape[1] or len(train)-X.shape[1]<3:return None
 return X,T,sd
models=[];pred=[];coverage=[]
scopes=[('discovery23',v[v.dataset.ne('Ogden2025')]),('all38',v)]+[(k,u) for k,u in v.groupby('dataset')]
for outcome in S['modifiers']:
 cols=[c for c in S['base_covariates'] if c!=outcome]
 for scope,u in scopes:
  u=u.reset_index(drop=True)
  for p in predictors:
   mat=matrices(u,u,cols,p)
   if mat is not None and len(u)-mat[0].shape[1]>=5:
    X,T,sd=mat;fit=sm.OLS(u[outcome],X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc[p]
    models.append(dict(scope=scope,outcome=outcome,predictor=p,n=len(u),df=fit.df_resid,beta=fit.params[p]/u[outcome].std(),slope=fit.params[p]/sd[p],CI_low=ci.iloc[0]/sd[p],CI_high=ci.iloc[1]/sd[p],P=fit.pvalues[p],family_target='four_states_'+scope,family_benchmark='benchmark27_'+scope+'_'+outcome))
   count=0
   for i in range(len(u)):
    tr=u.drop(index=i);te=u.iloc[[i]]
    b=matrices(tr,te,cols);a=matrices(tr,te,cols,p)
    if b is None or a is None:continue
    bp=b[1].to_numpy()@np.linalg.lstsq(b[0],tr[outcome],rcond=None)[0]
    ap=a[1].to_numpy()@np.linalg.lstsq(a[0],tr[outcome],rcond=None)[0]
    pred.append(dict(validation='leave_one_patient',scope=scope,outcome=outcome,predictor=p,patient=te.patient.iloc[0],observed=te[outcome].iloc[0],baseline=bp[0],augmented=ap[0]));count+=1
   coverage.append(dict(validation='leave_one_patient',scope=scope,outcome=outcome,predictor=p,expected=len(u),observed=count))
 for target in sorted(v.dataset.unique()):
  tr=v[v.dataset.ne(target)].reset_index(drop=True);te=v[v.dataset.eq(target)].reset_index(drop=True)
  for p in predictors:
   b=matrices(tr,te,cols,study=False);a=matrices(tr,te,cols,p,study=False)
   if b is None or a is None:
    coverage.append(dict(validation='leave_one_cohort',scope=target,outcome=outcome,predictor=p,expected=len(te),observed=0));continue
   bp=b[1].to_numpy()@np.linalg.lstsq(b[0],tr[outcome],rcond=None)[0];ap=a[1].to_numpy()@np.linalg.lstsq(a[0],tr[outcome],rcond=None)[0]
   for i in range(len(te)):pred.append(dict(validation='leave_one_cohort',scope=target,outcome=outcome,predictor=p,patient=te.patient.iloc[i],observed=te[outcome].iloc[i],baseline=bp[i],augmented=ap[i]))
   coverage.append(dict(validation='leave_one_cohort',scope=target,outcome=outcome,predictor=p,expected=len(te),observed=len(te)))
 print('STATE_DONE',outcome,flush=True)
t=pd.DataFrame(models);t['FDR_target']=np.nan;t['FDR_benchmark']=np.nan
for family,ids in t[t.predictor.eq('CGNL1')].groupby('family_target').groups.items():
 ids=list(ids);p=t.loc[ids,'P'].to_numpy();t.loc[ids,'FDR_target']=multipletests(np.r_[p,np.ones(4-len(p))],method='fdr_bh')[1][:len(p)]
for family,ids in t.groupby('family_benchmark').groups.items():
 ids=list(ids);p=t.loc[ids,'P'].to_numpy();t.loc[ids,'FDR_benchmark']=multipletests(np.r_[p,np.ones(27-len(p))],method='fdr_bh')[1][:len(p)]
t.to_csv(O/'E5_all_state_models.tsv',sep='\t',index=False)
t[t.predictor.eq('CGNL1')].to_csv(O/'E5_CGNL1_state_models.tsv',sep='\t',index=False)
pr=pd.DataFrame(pred);pr.to_csv(O/'E5_state_predictions.tsv',sep='\t',index=False)
perf=[]
for key,g in pr.groupby(['validation','scope','outcome','predictor']):
 b=(g.observed-g.baseline)**2;a=(g.observed-g.augmented)**2
 perf.append(dict(zip(['validation','scope','outcome','predictor'],key),n=len(g),MSE_baseline=b.mean(),MSE_augmented=a.mean(),improvement=1-a.mean()/b.mean(),patients_improved=int((a<b).sum())))
pf=pd.DataFrame(perf)
pf['rank']=pf.groupby(['validation','scope','outcome']).improvement.rank(method='min',ascending=False)
pf.to_csv(O/'E5_state_predictive_increment.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'E5_prediction_coverage.tsv',sep='\t',index=False)
old=pd.read_csv(O/'E2_CGNL1_state_links.tsv',sep='\t')
new=t[t.predictor.eq('CGNL1')]
errs=[]
for _,row in old.iterrows():
 q=new[new.scope.eq(row.dataset)&new.outcome.eq(row.state)].iloc[0]
 errs.extend([abs(q.beta-row.beta),abs(q.P-row.P)])
assert max(errs)<1e-9
(R/'E5_QA.json').write_text(json.dumps(dict(max_prior_reproduction_error=max(errs),all_prediction_folds_complete=all(x['expected']==x['observed'] for x in coverage),all_pass=True),indent=2))
print('STATE_INCREMENT_COMPLETE',flush=True)
