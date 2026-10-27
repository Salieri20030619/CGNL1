from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];B=R.parent;L=B/'threshold_ladder_20261002_v1';D=R/'results';C=json.loads((R/'SPECIFICATION.json').read_text());rng=np.random.default_rng(20261002)
# Quantify measurement and identities without adding exclusions.
rows=[];pieces=[]
for rec in json.loads((L/'CACHE_INDEX.json').read_text()):
 if rec['adjusted'] or not rec['infer']:continue
 with np.load(rec['path']) as f:z={k:f[k] for k in f.files if k!='barcode'}
 for stage in C['stages']:
  s=next(x for x in json.loads((L/'SPECIFICATION.json').read_text())['raw_stages'] if x['stage']==stage);mask=(z['qlib']>=s['umi'])&(z['qng']>=s['genes'])&(z['qmt']<=s['mt'])&(z['hit']>=s['hits'])&(z['margin']>=s['margin']);n=int(mask.sum())
  if n<s['mincells']:continue
  raw=z['counts'][:,mask];lib=z['lib'][mask];row=dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],stage=stage,cells=n,CGNL1_total=float(raw[0].sum()),CGNL1_positive_cells=int((raw[0]>0).sum()),CGNL1_detection_fraction=float((raw[0]>0).mean()),CGNL1_CP10K=float(raw[0].sum()/lib.sum()*1e4),median_UMI=float(np.median(lib)),median_genes=float(np.median(z['qng'][mask])))
  for j,g in enumerate(['CGNL1','CLDN5','ESAM','TJP1','OCLN']):row[g+'_counts']=float(raw[j].sum())
  rows.append(row)
 q=(z['qlib']>=500)&(z['qng']>=200)&(z['qmt']<=25)
 for group,k in [('one_marker_removed',q&(z['hit']==1)),('two_plus_retained',q&(z['hit']>=2))]:
  n=int(k.sum())
  if not n:continue
  c=z['counts'][:,k];lib=z['lib'][k];val=np.log1p(c/lib*1e4);pieces.append(dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],group=group,cells=n,median_UMI=float(np.median(lib)),median_genes=float(np.median(z['ng'][k])),median_identity_margin=float(np.median(z['margin'][k])),CGNL1_positive_fraction=float((c[0]>0).mean()),CGNL1_total=float(c[0].sum()),CGNL1_mean_logCP10K=float(val[0].mean()),junction_mean_logCP10K=float(val[1:].mean())))
rows=pd.DataFrame(rows);rows.to_csv(D/'measurement_by_lesion.tsv',sep='\t',index=False);pr=rows.groupby(['stage','dataset','donor'],as_index=False).mean(numeric_only=True);pr.to_csv(D/'measurement_by_patient.tsv',sep='\t',index=False)
agg=pr.groupby(['stage','dataset'],as_index=False).agg(n_patients=('donor','nunique'),zero_CGNL1_patients=('CGNL1_total',lambda x:int((x==0).sum())),median_CGNL1_counts=('CGNL1_total','median'),median_CGNL1_detection=('CGNL1_detection_fraction','median'),median_UMI=('median_UMI','median'),median_endo_cells=('cells','median'));agg.to_csv(D/'measurement_by_cohort.tsv',sep='\t',index=False)
pieces=pd.DataFrame(pieces);pieces.to_csv(D/'removed_retained_cell_profiles.tsv',sep='\t',index=False)
pp=pieces.groupby(['dataset','donor','group'],as_index=False).mean(numeric_only=True);pp.to_csv(D/'removed_retained_patient_profiles.tsv',sep='\t',index=False)
p=pd.read_csv(L/'results/patient_profiles_all_stages.tsv',sep='\t');p=p[(p.branch=='raw')&p.infer];models=[];influence=[];heter=[]
def design(v):
 Z=pd.get_dummies(v.dataset,drop_first=True,dtype=float);Z=sm.add_constant(Z,has_constant='add');F=Z.copy();F['CGNL1']=v.CGNL1;return Z.to_numpy(),F
for stage in C['primary_diagnostic_stages']:
 for scope in ['all','external']:
  v=p[p.stage==stage].copy()
  if scope=='external':v=v[v.dataset!='GSE225857']
  v=v.reset_index(drop=True);Z,F=design(v);y=v.junction.to_numpy();x=v.CGNL1.to_numpy();assert np.linalg.matrix_rank(F.to_numpy())==F.shape[1]
  ols=sm.OLS(y,F).fit();hc=sm.OLS(y,F).fit(cov_type='HC3',use_t=True)
  xr=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];yr=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0];stat=float(xr@yr);extreme=0;groups=[np.flatnonzero(v.dataset.to_numpy()==ds) for ds in v.dataset.unique()]
  for i in range(4999):
   yp=yr.copy()
   for ids in groups:yp[ids]=rng.permutation(yr[ids])
   extreme += abs(float(xr@yp))>=abs(stat)-1e-12
  perm=(extreme+1)/5000
  for method,mod,pval in [('OLS',ols,float(ols.pvalues['CGNL1'])),('HC3',hc,float(hc.pvalues['CGNL1'])),('within_study_permutation',ols,perm)]:
   ci=mod.conf_int().loc['CGNL1'];models.append(dict(stage=stage,scope=scope,method=method,n=len(v),beta=float(mod.params['CGNL1'])*x.std(ddof=1)/y.std(ddof=1),slope=float(mod.params['CGNL1']),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),CI_method='OLS' if method=='within_study_permutation' else method,P=pval,family='method_robustness_'+scope))
  # Study interactions diagnose heterogeneous slopes, not a reason to discard a study.
  full=F.copy()
  for ds in v.dataset.unique()[1:]:full['CGNL1_x_'+ds]=v.CGNL1*(v.dataset==ds)
  if np.linalg.matrix_rank(full.to_numpy())==full.shape[1] and len(v)-full.shape[1]>=3:
   fm=sm.OLS(y,full).fit();fv,pv,df=fm.compare_f_test(ols);heter.append(dict(stage=stage,scope=scope,n=len(v),F=float(fv),P=float(pv),df_difference=float(df),residual_df=float(fm.df_resid),family='study_slope_heterogeneity_'+scope))
  if scope=='all':
   for ds in v.dataset.unique():
    w=v[v.dataset!=ds].copy();_,W=design(w);mod=sm.OLS(w.junction,W).fit();ci=mod.conf_int().loc['CGNL1'];influence.append(dict(stage=stage,omitted_study=ds,n=len(w),beta=float(mod.params['CGNL1'])*w.CGNL1.std()/w.junction.std(),slope=float(mod.params['CGNL1']),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(mod.pvalues['CGNL1']),family='leave_one_study_'+stage,scope='diagnostic_not_cohort_selection'))
for name,data in [('method_robustness',models),('study_slope_heterogeneity',heter),('leave_one_study_out',influence)]:
 t=pd.DataFrame(data);t['FDR']=np.nan
 for fam,ids in t.groupby('family').groups.items():t.loc[ids,'FDR']=multipletests(t.loc[ids,'P'],method='fdr_bh')[1]
 t['nominal_positive_association']=(t.P<.05)&(t.beta>0) if 'beta' in t else t.P<.05;t.to_csv(D/(name+'.tsv'),sep='\t',index=False)
summary={'exit_code':0,'finished':datetime.datetime.now().astimezone().isoformat(),'formal_positive_gate':'raw target-effect P<0.05; FDR retained','method_robustness_tests':len(models),'study_interaction_tests':len(heter),'study_influence_models':len(influence),'removed_retained_patient_profiles':len(pp),'no_original_data_modified':True}
(R/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
