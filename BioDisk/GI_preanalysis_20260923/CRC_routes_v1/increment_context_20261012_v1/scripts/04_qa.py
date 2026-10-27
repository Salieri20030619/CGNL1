from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib
import numpy as np,pandas as pd,statsmodels.formula.api as smf
from scipy import stats
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'increment_context_20261012_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
V=pd.read_csv(O/'patient_profiles_all_variants.tsv',sep='\t')
checks={};checks['frozen_spec']=hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
checks['postprimary_amendment_frozen']=hashlib.sha256((R/'POST_PRIMARY_STATE_INCREMENT.json').read_bytes()).hexdigest()==json.loads((R/'POST_PRIMARY_STATE_FREEZE.json').read_text())['sha256']
errors=[]
for _,row in pd.read_csv(O/'E1_measurement_slopes.tsv',sep='\t').iterrows():
 variant='legacy' if row.variant=='within_cohort_rank_legacy' else row.variant
 v=V[V.variant.eq(variant)&V.definition.isin(['discovery_frozen',row.definition])].copy()
 if row.variant=='within_cohort_rank_legacy':
  for col in ['CGNL1','actomyosin']+S['base_covariates']:
   v[col]=v.groupby('dataset')[col].transform(lambda x:stats.rankdata(x)/(len(x)+1))
 if row.scope=='discovery23':v=v[v.dataset.ne('Ogden2025')]
 else:v=v[v.dataset.eq(row.scope)]
 assert v.patient.is_unique
 fit=smf.ols('actomyosin ~ CGNL1 + depth + epithelial_identity + mature_colonocyte + C(dataset)',data=v).fit(cov_type='HC3',use_t=True)
 errors.extend([abs(fit.params.CGNL1*v.CGNL1.std()/v.actomyosin.std()-row.beta),abs(fit.pvalues.CGNL1-row.P)])
checks['independent_formula_max_error']=max(errors)
assert max(errors)<1e-8
def bh(p,m):
 p=np.array(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*m/np.arange(1,len(p)+1))[::-1])[::-1]);return q
fdr_err=[]
for file,group,p,q,m in [
 ('E1_measurement_slopes.tsv','family','P','FDR',4),('E1_heterogeneity.tsv','family','P','FDR',2),
 ('E2_state_interactions.tsv','family','P','FDR',4),('E2_CGNL1_state_links.tsv','family','P','FDR',4),
 ('E4_benchmark_models.tsv','family','P','FDR',27),('E5_CGNL1_state_models.tsv','family_target','P','FDR_target',4),
 ('E5_all_state_models.tsv','family_benchmark','P','FDR_benchmark',27)]:
 t=pd.read_csv(O/file,sep='\t')
 for _,g in t.groupby(group):
  good=g[p].notna();g=g[good]
  fdr_err.extend(abs(bh(g[p],m)-g[q]).tolist())
checks['independent_BH_max_error']=max(fdr_err);assert max(fdr_err)<1e-10
# Verify full scores directly from gene-level columns on every measurement row.
errs=[]
for name,genes in {**S['panels'],**S['outcomes']}.items():
 if all(g in V for g in genes):errs.extend(abs(V[genes].mean(1)-V[name]).dropna().tolist())
checks['score_reconstruction_error']=max(errs);assert max(errs)<1e-10
# Recompute all reported errors from preserved held-out predictions.
for file,summary,keys,aug,imp in [
 ('E3_heldout_predictions.tsv','E3_prediction_performance.tsv',['variant','validation','scope','modifier'],'additive','additive_improvement'),
 ('E5_state_predictions.tsv','E5_state_predictive_increment.tsv',['validation','scope','outcome','predictor'],'augmented','improvement')]:
 p=pd.read_csv(O/file,sep='\t');p['modifier']=p.get('modifier',pd.Series('none',index=p.index)).fillna('none')
 t=pd.read_csv(O/summary,sep='\t');t['modifier']=t.get('modifier',pd.Series('none',index=t.index)).fillna('none')
 assert not p.duplicated(keys+['patient']).any()
 errs=[]
 for key,g in p.groupby(keys):
  match=np.ones(len(t),dtype=bool)
  for k,val in zip(keys,key):match &= t[k].eq(val).to_numpy()
  row=t.loc[match].iloc[0]
  computed=1-((g[aug]-g.observed)**2).mean()/((g.baseline-g.observed)**2).mean()
  errs.append(abs(computed-row[imp]))
 checks[file+'_max_metric_error']=max(errs);assert max(errs)<1e-8
# Explicitly enumerate absent prediction branches; never interpret them as zero effects.
t=pd.read_csv(O/'E3_prediction_performance.tsv',sep='\t').fillna({'modifier':'none'})
coverage=[]
nmap={'GSE245552':13,'GSE285991':10,'Ogden2025':15,'discovery23':23,'all38':38}
for variant in ['legacy','pooled_common_genes']:
 for validation,scopes in [('leave_one_patient',list(nmap)),('leave_one_cohort',['GSE245552','GSE285991','Ogden2025'])]:
  for scope in scopes:
   for modifier in ['none']+S['modifiers']:
    g=t[t.variant.eq(variant)&t.validation.eq(validation)&t.scope.eq(scope)&t.modifier.eq(modifier)]
    n=int(g.n.iloc[0]) if len(g) else 0
    coverage.append(dict(variant=variant,validation=validation,scope=scope,modifier=modifier,expected=nmap[scope],observed=n,status='complete' if n==nmap[scope] else 'not_estimable_df_rank' if n==0 else 'partial'))
pd.DataFrame(coverage).to_csv(O/'E3_prediction_coverage.tsv',sep='\t',index=False)
checks['E3_unestimable_branches']=sum(x['observed']==0 for x in coverage)
checks['E3_partial_branches']=sum(0<x['observed']<x['expected'] for x in coverage)
assert checks['E3_partial_branches']==0
# Feature-only standardization is equivalent to raw design OLS. Independently
# reconstruct every CGNL1 four-state LOPO/LOCO prediction with raw formula design.
import patsy
v=V[V.variant.eq('pooled_common_genes')&V.definition.isin(['discovery_frozen','author_epithelial_QC'])].reset_index(drop=True)
pr=pd.read_csv(O/'E5_state_predictions.tsv',sep='\t');pr=pr[pr.predictor.eq('CGNL1')]
errs=[]
for key,g in pr.groupby(['validation','scope','outcome']):
 validation,scope,outcome=key
 if validation=='leave_one_patient':
  u=v if scope=='all38' else v[v.dataset.ne('Ogden2025')] if scope=='discovery23' else v[v.dataset.eq(scope)]
 else:u=v
 cols=[c for c in S['base_covariates'] if c!=outcome]
 formula=outcome+' ~ '+' + '.join(cols+['CGNL1']+(['C(dataset)'] if validation=='leave_one_patient' else []))
 for _,row in g.iterrows():
  tr=u[u.patient.ne(row.patient)] if validation=='leave_one_patient' else u[u.dataset.ne(scope)]
  te=u[u.patient.eq(row.patient)]
  fit=smf.ols(formula,data=tr).fit()
  errs.append(abs(float(fit.predict(te).iloc[0])-row.augmented))
checks['independent_state_prediction_max_error']=max(errs);assert max(errs)<1e-8
# Source fingerprints establish that historical artifacts were not modified.
source_ok=[]
for item in json.loads((R/'SOURCE_FINGERPRINTS.json').read_text()):
 p=Path(item['path']);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 source_ok.append(h.hexdigest()==item['sha256'])
checks['all_prior_sources_unchanged']=all(source_ok);assert all(source_ok)
checks['all_pass']=True
(R/'NUMERICAL_QA.json').write_text(json.dumps(checks,indent=2))
print('QA_COMPLETE',json.dumps(checks),flush=True)
