from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,hashlib,datetime
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'new_cohorts_20261013_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
v=pd.read_csv(O/'primary_patient_profiles.tsv',sep='\t')
eligible=[];models=[];pred=[];num=[];transfer=[];bench=[]
def matrix(train,test,cols,predictor=None):
 features=list(cols)+([predictor] if predictor else [])
 mu=train[features].mean();sd=train[features].std();keep=sd[sd>1e-10].index
 if predictor and predictor not in keep:return None
 X=pd.DataFrame({'const':np.ones(len(train))},index=train.index);T=pd.DataFrame({'const':np.ones(len(test))},index=test.index)
 for c in keep:X[c]=(train[c]-mu[c])/sd[c];T[c]=(test[c]-mu[c])/sd[c]
 if not np.isfinite(X.to_numpy()).all() or not np.isfinite(T.to_numpy()).all() or np.linalg.matrix_rank(X)<X.shape[1] or np.linalg.cond(X)>1e8:return None
 return X,T,sd
for (dataset,definition),g in v.groupby(['dataset','definition']):
 for minimum in [20,50]:
  u=g[g.cells.ge(minimum)].reset_index(drop=True)
  for model,cols in S['models'].items():
   for outcome in S['outcomes']:
    z=u.dropna(subset=cols+[outcome,'CGNL1']).reset_index(drop=True)
    mat=matrix(z,z,cols,'CGNL1') if len(z)>1 else None
    qualified=len(z)>=10 and mat is not None and len(z)-mat[0].shape[1]>=5 and z[outcome].std()>1e-10
    info=dict(dataset=dataset,definition=definition,min_cells=minimum,model=model,outcome=outcome,n=len(z))
    eligible.append(dict(info,status='estimated' if qualified else 'not_estimable_sample_rank_df'))
    if not qualified:continue
    X,T,sd=mat;fit=sm.OLS(z[outcome],X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['CGNL1']/sd.CGNL1
    models.append(dict(info,beta=fit.params.CGNL1/z[outcome].std(),slope=fit.params.CGNL1/sd.CGNL1,CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=fit.pvalues.CGNL1,df=fit.df_resid,family=dataset+'_'+definition+'_'+str(minimum)+'_'+model+'_two_endpoints'))
    N=X.drop(columns='CGNL1').to_numpy();xx=X.CGNL1.to_numpy();yy=z[outcome].to_numpy();rx=xx-N@np.linalg.lstsq(N,xx,rcond=None)[0];ry=yy-N@np.linalg.lstsq(N,yy,rcond=None)[0]
    num.append(abs(float(rx@ry/(rx@rx))-fit.params.CGNL1))
    for i in range(len(z)):
     tr=z.drop(index=i);te=z.iloc[[i]];b=matrix(tr,te,cols);a=matrix(tr,te,cols,'CGNL1')
     if b is None or a is None or len(tr)-a[0].shape[1]<3:continue
     bp=float((b[1].to_numpy()@np.linalg.lstsq(b[0],tr[outcome],rcond=None)[0])[0]);ap=float((a[1].to_numpy()@np.linalg.lstsq(a[0],tr[outcome],rcond=None)[0])[0])
     pred.append(dict(info,patient=te.patient.iloc[0],observed=te[outcome].iloc[0],baseline=bp,augmented=ap))
 # Fixed six ordinary markers, same main20 M1 baseline.
 if definition=='main':
  z=g[g.cells.ge(20)].reset_index(drop=True)
  for outcome in S['outcomes']:
   for predictor in ['CGNL1']+S['ordinary_markers']:
    cols=S['models']['M1'];mat=matrix(z,z,cols,predictor)
    if len(z)<10 or mat is None or len(z)-mat[0].shape[1]<5:continue
    X,_,sd=mat;fit=sm.OLS(z[outcome],X).fit(cov_type='HC3',use_t=True)
    for i in range(len(z)):
     tr=z.drop(index=i);te=z.iloc[[i]];b=matrix(tr,te,cols);a=matrix(tr,te,cols,predictor)
     if b is None or a is None:continue
     bp=float((b[1].to_numpy()@np.linalg.lstsq(b[0],tr[outcome],rcond=None)[0])[0]);ap=float((a[1].to_numpy()@np.linalg.lstsq(a[0],tr[outcome],rcond=None)[0])[0])
     bench.append(dict(dataset=dataset,outcome=outcome,predictor=predictor,patient=te.patient.iloc[0],observed=te[outcome].iloc[0],baseline=bp,augmented=ap))
t=pd.DataFrame(models);t['FDR']=np.nan
for _,ids in t.groupby('family').groups.items():
 ids=list(ids);p=t.loc[ids,'P'].to_numpy();t.loc[ids,'FDR']=multipletests(np.r_[p,np.ones(2-len(p))],method='fdr_bh')[1][:len(p)]
t.to_csv(O/'N1_CGNL1_models.tsv',sep='\t',index=False);pd.DataFrame(eligible).to_csv(O/'N1_model_eligibility.tsv',sep='\t',index=False)
p=pd.DataFrame(pred);p.to_csv(O/'N1_heldout_predictions.tsv',sep='\t',index=False)
perf=[]
keys=['dataset','definition','min_cells','model','outcome']
for key,g in p.groupby(keys):
 b=(g.baseline-g.observed)**2;a=(g.augmented-g.observed)**2
 perf.append(dict(zip(keys,key),n=len(g),MSE_baseline=b.mean(),MSE_augmented=a.mean(),improvement=1-a.mean()/b.mean(),patients_improved=int((a<b).sum())))
pd.DataFrame(perf).to_csv(O/'N1_predictive_increment.tsv',sep='\t',index=False)
bp=pd.DataFrame(bench);bp.to_csv(O/'N2_benchmark_predictions.tsv',sep='\t',index=False);br=[]
for key,g in bp.groupby(['dataset','outcome','predictor']):
 b=(g.baseline-g.observed)**2;a=(g.augmented-g.observed)**2
 br.append(dict(zip(['dataset','outcome','predictor'],key),n=len(g),improvement=1-a.mean()/b.mean()))
pd.DataFrame(br).to_csv(O/'N2_benchmark_increment.tsv',sep='\t',index=False)
# Fixed old-discovery transfer and new-primary-cohort transfer; no target fitting.
old=pd.read_csv(B/'increment_context_20261012_v1/results/patient_profiles_all_variants.tsv',sep='\t')
old=old[old.variant.eq('pooled_allgenes')&old.definition.eq('discovery_frozen')].reset_index(drop=True)
cohorts={'old_discovery23':old}
cohorts.update({k:g.reset_index(drop=True) for k,g in v[v.definition.eq('main')].groupby('dataset')})
for train_name,train in cohorts.items():
 for test_name,test in cohorts.items():
  if train_name==test_name or test_name=='old_discovery23':continue
  for outcome in S['outcomes']:
   cols=S['models']['M1'];b=matrix(train,test,cols);a=matrix(train,test,cols,'CGNL1')
   if len(train)<10 or b is None or a is None or len(train)-a[0].shape[1]<5:continue
   bpred=b[1].to_numpy()@np.linalg.lstsq(b[0],train[outcome],rcond=None)[0];apred=a[1].to_numpy()@np.linalg.lstsq(a[0],train[outcome],rcond=None)[0]
   for i in range(len(test)):transfer.append(dict(train=train_name,test=test_name,outcome=outcome,patient=test.patient.iloc[i],observed=test[outcome].iloc[i],baseline=bpred[i],augmented=apred[i]))
p=pd.DataFrame(transfer);p.to_csv(O/'N3_transfer_predictions.tsv',sep='\t',index=False)
tr=[]
for key,g in p.groupby(['train','test','outcome']):
 b=(g.baseline-g.observed)**2;a=(g.augmented-g.observed)**2
 tr.append(dict(zip(['train','test','outcome'],key),n=len(g),MSE_baseline=b.mean(),MSE_augmented=a.mean(),improvement=1-a.mean()/b.mean(),patients_improved=int((a<b).sum())))
pd.DataFrame(tr).to_csv(O/'N3_transfer_performance.tsv',sep='\t',index=False)
# Small liver cohorts: preserve all patient effects, no multivariable inference.
l=pd.read_csv(O/'liver_patient_profiles.tsv',sep='\t');desc=[]
for (ds,site),g in l.groupby(['dataset','site']):
 for outcome in S['outcomes']:
  desc.append(dict(dataset=ds,site=site,outcome=outcome,n=len(g),spearman_r=float(stats.spearmanr(g.CGNL1,g[outcome]).statistic) if len(g)>=3 else np.nan,P=np.nan,FDR=np.nan,evidence='small-cohort descriptive, no adjusted incremental validation'))
pd.DataFrame(desc).to_csv(O/'N4_liver_descriptive.tsv',sep='\t',index=False)
pairs=[]
for ds,g in l.groupby('dataset'):
 if ds!='GSE178318':continue
 for patient,u in g.groupby('patient'):
  if set(u.site)!={'primary_CRC','liver_metastasis'}:continue
  a=u[u.site.eq('primary_CRC')].iloc[0];b=u[u.site.eq('liver_metastasis')].iloc[0]
  for feature in ['CGNL1','actomyosin','cycling']:pairs.append(dict(dataset=ds,patient=patient,feature=feature,primary=a[feature],liver=b[feature],delta=b[feature]-a[feature]))
pd.DataFrame(pairs).to_csv(O/'N4_paired_site_descriptive.tsv',sep='\t',index=False)
qa=dict(max_FWL_error=max(num),models=len(t),prediction_rows=len(pred),all_pass=bool(max(num)<1e-8))
assert qa['all_pass'];(R/'ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2))
print('ANALYSIS_COMPLETE',json.dumps(qa),flush=True)
print(t[t.definition.eq('main')&t.min_cells.eq(20)&t.model.eq('M1')].to_string(index=False))
print(pd.DataFrame(perf).query("definition=='main' and min_cells==20 and model=='M1'").to_string(index=False))
