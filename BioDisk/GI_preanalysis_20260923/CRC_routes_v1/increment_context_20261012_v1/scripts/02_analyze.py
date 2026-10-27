from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib,datetime
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'increment_context_20261012_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
V=pd.read_csv(O/'patient_profiles_all_variants.tsv',sep='\t')
base=S['base_covariates'];mods=S['modifiers']
controls=json.loads((B/'epithelial_independence_20261012_v1/results/A_CONTROL_FREEZE.json').read_text())['controls']
genes=pd.read_csv(B/'epithelial_independence_20261012_v1/results/A_patient_allgenes_logCP10k.tsv.gz',sep='\t',index_col=0)
for p in controls:
 ix=V.variant.eq('legacy')&V.dataset.ne('Ogden2025')
 V.loc[ix,p]=V.loc[ix,'patient'].map(genes[p])
elig=[];fwl=[]
def get(variant='legacy',definition='author_epithelial_QC'):
 d=V[V.variant.eq('legacy' if variant=='within_cohort_rank_legacy' else variant)&V.definition.isin(['discovery_frozen',definition])].copy().reset_index(drop=True)
 assert d.patient.is_unique and len(d)==(38 if definition=='author_epithelial_QC' else 36)
 if variant=='within_cohort_rank_legacy':
  for col in ['CGNL1','actomyosin']+list(dict.fromkeys(base+mods)):
   d[col]=d.groupby('dataset')[col].transform(lambda x:stats.rankdata(x)/(len(x)+1))
 return d
def design(v,cols,target='CGNL1',study=True):
 d=pd.DataFrame({'const':np.ones(len(v))},index=v.index)
 for col in cols:d[col]=v[col].astype(float)
 if study:
  for k in sorted(v.dataset.unique())[1:]:d['study_'+k]=v.dataset.eq(k).astype(float)
 if target:d['target']=v[target].astype(float)
 return d
def estimate(v,D,outcome='actomyosin',term='target',label=''):
 y=v[outcome].astype(float)
 if not np.isfinite(D.to_numpy()).all() or not np.isfinite(y).all():
  elig.append(dict(label=label,n=len(v),status='missing'));return None
 active=[c for c in D if c=='const' or D[c].std()>1e-10]
 D=D[active].copy()
 # Scale columns for stable algebra, retaining raw-scale coefficients in output.
 mu=D.mean();sd=D.std();mu.loc['const']=0;sd.loc['const']=1
 D=(D-mu)/sd
 df=len(v)-D.shape[1]
 if term not in D or df<5 or np.linalg.matrix_rank(D)<D.shape[1] or y.std()<1e-10 or np.linalg.cond(D)>1e8:
  elig.append(dict(label=label,n=len(v),status='rank_df_variance',df=df));return None
 fit=sm.OLS(y,D).fit(cov_type='HC3',use_t=True)
 ci=fit.conf_int().loc[term]/sd[term]
 slope=float(fit.params[term]/sd[term])
 N=D.drop(columns=term).to_numpy();x=D[term].to_numpy()
 rx=x-N@np.linalg.lstsq(N,x,rcond=None)[0];ry=y-N@np.linalg.lstsq(N,y,rcond=None)[0]
 err=abs(float(rx@ry/(rx@rx))-float(fit.params[term]))
 fwl.append(err)
 row=dict(n=len(v),df=df,slope=slope,beta=float(fit.params[term]/y.std()),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(fit.pvalues[term]),max_leverage=float(sm.OLS(y,D).fit().get_influence().hat_matrix_diag.max()))
 elig.append(dict(label=label,n=len(v),status='estimated',df=df))
 return row,fit,D,sd
def adj(t,groups,m=None):
 t=t.copy();t['FDR']=np.nan
 if not len(t):return t
 for _,ids in t.groupby(groups,dropna=False).groups.items():
  ids=list(ids);ok=[i for i in ids if np.isfinite(t.loc[i,'P'])]
  if ok:
   pp=t.loc[ok,'P'].to_numpy();M=len(ids) if m is None else m
   # Missing hypotheses remain in their prespecified family.
   t.loc[ok,'FDR']=multipletests(np.r_[pp,np.ones(max(0,M-len(pp)))],method='fdr_bh')[1][:len(pp)]
 return t
def scopes(d):
 return [('discovery23',d[d.dataset.ne('Ogden2025')].reset_index(drop=True))]+[(k,v.reset_index(drop=True)) for k,v in d.groupby('dataset')]
rows=[];hets=[];residuals=[]
for definition in ['author_epithelial_QC','frozen_identity_intersection']:
 for variant in S['measurement_variants']:
  d=get(variant,definition)
  for scope,v in scopes(d):
   est=estimate(v,design(v,base),label='measurement:'+definition+':'+variant+':'+scope)
   if est:rows.append(dict(definition=definition,variant=variant,scope=scope,family='measurement_slopes_'+definition+'_'+variant,**est[0]))
  # A fully cohort-specific nuisance model compares standardized target slopes.
  D=pd.DataFrame(index=d.index);datasets=sorted(d.dataset.unique())
  zy=d.groupby('dataset').actomyosin.transform(lambda x:(x-x.mean())/x.std())
  for k in datasets:
   flag=d.dataset.eq(k).astype(float);D['intercept_'+k]=flag
   for col in base+['CGNL1']:
    z=(d[col]-d.loc[d.dataset.eq(k),col].mean())/d.loc[d.dataset.eq(k),col].std()
    D[col+'_'+k]=flag*z
  if np.linalg.matrix_rank(D)==D.shape[1] and len(d)-D.shape[1]>=5:
   fit=sm.OLS(zy,D).fit(cov_type='HC3',use_t=True)
   contrasts=[]
   C=np.zeros((2,D.shape[1]))
   for j in range(2):
    C[j,D.columns.get_loc('CGNL1_'+datasets[j])]=1;C[j,D.columns.get_loc('CGNL1_'+datasets[2])]=-1
   contrasts.append(('all_three_slopes',C))
   c=np.zeros((1,D.shape[1]))
   for k in datasets:c[0,D.columns.get_loc('CGNL1_'+k)]=1 if k=='Ogden2025' else -.5
   contrasts.append(('external_vs_equal_discovery_mean',c))
   for name,C in contrasts:
    test=fit.f_test(C)
    hets.append(dict(definition=definition,variant=variant,test=name,n=len(d),F=float(test.fvalue),df_num=float(test.df_num),df_den=float(test.df_denom),P=float(test.pvalue),contrast=float((C@fit.params)[0]) if len(C)==1 else np.nan,family='heterogeneity_'+definition+'_'+variant))
  else:elig.append(dict(label='heterogeneity:'+definition+':'+variant,n=len(d),status='not_estimable'))
  if variant=='pooled_common_genes' and definition=='author_epithelial_QC':
   for scope,v in scopes(d):
    N=design(v,base,target=None).to_numpy()
    rx=v.CGNL1-N@np.linalg.lstsq(N,v.CGNL1,rcond=None)[0];ry=v.actomyosin-N@np.linalg.lstsq(N,v.actomyosin,rcond=None)[0]
    for i in range(len(v)):residuals.append(dict(scope=scope,patient=v.patient.iloc[i],dataset=v.dataset.iloc[i],x=float(rx.iloc[i]),y=float(ry.iloc[i])))
a=adj(pd.DataFrame(rows),['family'],4);a.to_csv(O/'E1_measurement_slopes.tsv',sep='\t',index=False)
adj(pd.DataFrame(hets),['family'],2).to_csv(O/'E1_heterogeneity.tsv',sep='\t',index=False)
pd.DataFrame(residuals).to_csv(O/'E1_residuals.tsv',sep='\t',index=False)
print('E1_COMPLETE',flush=True)
# State interaction: cohort-specific standardized nuisance and target effects.
d=get('pooled_common_genes')
inter=[];links=[];state_desc=[]
for k,v in d.groupby('dataset'):
 for col in ['CGNL1','actomyosin']+base+mods:
  state_desc.append(dict(dataset=k,variable=col,n=len(v),mean=v[col].mean(),sd=v[col].std(),minimum=v[col].min(),maximum=v[col].max()))
for modifier in mods:
 z=d.copy()
 for col in list(dict.fromkeys(base+['CGNL1','actomyosin',modifier])):
  z[col]=z.groupby('dataset')[col].transform(lambda x:(x-x.mean())/x.std())
 D=pd.DataFrame(index=z.index)
 for k in sorted(z.dataset.unique()):
  f=z.dataset.eq(k).astype(float);D['intercept_'+k]=f
  for col in list(dict.fromkeys(base+[modifier,'CGNL1'])):D[col+'_'+k]=f*z[col]
 D['interaction']=z.CGNL1*z[modifier]
 # Direct fit retains all cohort intercepts; avoid extra intercept.
 df=len(z)-D.shape[1]
 if np.linalg.matrix_rank(D)==D.shape[1] and df>=5:
  fit=sm.OLS(z.actomyosin,D).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['interaction']
  inter.append(dict(modifier=modifier,n=len(z),df=df,beta=fit.params.interaction,CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=fit.pvalues.interaction,family='state_interactions_4'))
 else:inter.append(dict(modifier=modifier,n=len(z),df=df,beta=np.nan,P=np.nan,family='state_interactions_4'))
 for k,v in d.groupby('dataset'):
  v=v.reset_index(drop=True);cols=[x for x in base if x!=modifier]
  est=estimate(v,design(v,cols),outcome=modifier,label='state_link:'+k+':'+modifier)
  if est:links.append(dict(dataset=k,state=modifier,family='CGNL1_state_links_'+k,**est[0]))
adj(pd.DataFrame(inter),['family'],4).to_csv(O/'E2_state_interactions.tsv',sep='\t',index=False)
adj(pd.DataFrame(links),['family'],4).to_csv(O/'E2_CGNL1_state_links.tsv',sep='\t',index=False)
pd.DataFrame(state_desc).drop_duplicates().to_csv(O/'E2_distribution_summary.tsv',sep='\t',index=False)
# Therapy interaction uses15 patients, not cells or repeated subtype bins.
e=d[d.dataset.eq('Ogden2025')].reset_index(drop=True).copy();e['treated']=e.therapy.eq('Treated').astype(float)
D=design(e,base+['treated']);D['interaction']=(e.CGNL1-e.CGNL1.mean())*e.treated
est=estimate(e,D,term='interaction',label='therapy_interaction')
therapy=dict(naive_n=int(e.treated.eq(0).sum()),treated_n=int(e.treated.eq(1).sum()),status='estimated' if est else 'not_estimable',family='single_prespecified_therapy_interaction')
if est:therapy.update(est[0])
(R/'THERAPY_INTERACTION.json').write_text(json.dumps(therapy,indent=2))
tr=[]
for k,v in e.groupby('therapy'):
 tr.append(dict(therapy=k,n=len(v),unadjusted_slope=np.polyfit(v.CGNL1,v.actomyosin,1)[0],pearson_r=v.CGNL1.corr(v.actomyosin),P=np.nan,evidence='descriptive small subgroup; no multivariable subgroup inference'))
pd.DataFrame(tr).to_csv(O/'E2_therapy_descriptive.tsv',sep='\t',index=False)
print('E2_COMPLETE',flush=True)
# No held-out target or feature statistics enter model fitting.
def predict(train,test,cols,mode='baseline',modifier=None,study=True,cohort_targets=False,predictor='CGNL1'):
 cols=list(dict.fromkeys(cols+([modifier] if modifier else [])))
 num=list(dict.fromkeys(cols+([predictor] if mode!='baseline' else [])))
 mu=train[num].mean();sd=train[num].std();active=sd[sd.gt(1e-10)].index.tolist()
 if mode!='baseline' and predictor not in active:return None
 x=(train[active]-mu[active])/sd[active];t=(test[active]-mu[active])/sd[active]
 D=pd.DataFrame({'const':np.ones(len(train))},index=train.index);T=pd.DataFrame({'const':np.ones(len(test))},index=test.index)
 for col in active:D[col]=x[col];T[col]=t[col]
 if study:
  levels=sorted(train.dataset.unique())
  if not set(test.dataset)<=set(levels):return None
  for k in levels[1:]:
   D['study_'+k]=train.dataset.eq(k).astype(float);T['study_'+k]=test.dataset.eq(k).astype(float)
   if cohort_targets and mode!='baseline':
    D['target_'+k]=D['study_'+k]*x[predictor];T['target_'+k]=T['study_'+k]*t[predictor]
 if mode=='interaction':
  if modifier not in active:return None
  D['interaction']=x[predictor]*x[modifier];T['interaction']=t[predictor]*t[modifier]
 if not np.isfinite(D.to_numpy()).all() or not np.isfinite(T.to_numpy()).all() or len(train)-D.shape[1]<3 or np.linalg.matrix_rank(D)<D.shape[1] or np.linalg.cond(D)>1e8:return None
 return T.to_numpy()@np.linalg.lstsq(D,train.actomyosin,rcond=None)[0]
pred=[];overlap=[]
for variant in ['legacy','pooled_common_genes']:
 d=get(variant)
 for scope,v in scopes(d)+[('all38',d)]:
  v=v.reset_index(drop=True)
  for modifier in [None]+mods:
   for i in range(len(v)):
    train=v.drop(index=i);test=v.iloc[[i]]
    bp=predict(train,test,base,'baseline',modifier)
    ap=predict(train,test,base,'additive',modifier)
    ip=predict(train,test,base,'interaction',modifier) if modifier else ap
    if bp is None or ap is None or ip is None:continue
    pred.append(dict(variant=variant,validation='leave_one_patient',scope=scope,modifier=modifier or 'none',patient=test.patient.iloc[0],observed=test.actomyosin.iloc[0],baseline=bp[0],additive=ap[0],interaction=ip[0]))
 for target in sorted(d.dataset.unique()):
  tr=d[d.dataset.ne(target)].reset_index(drop=True);te=d[d.dataset.eq(target)].reset_index(drop=True)
  for modifier in [None]+mods:
   bp=predict(tr,te,base,'baseline',modifier,study=False);ap=predict(tr,te,base,'additive',modifier,study=False)
   ip=predict(tr,te,base,'interaction',modifier,study=False) if modifier else ap
   if bp is None or ap is None or ip is None:continue
   for i in range(len(te)):pred.append(dict(variant=variant,validation='leave_one_cohort',scope=target,modifier=modifier or 'none',patient=te.patient.iloc[i],observed=te.actomyosin.iloc[i],baseline=bp[i],additive=ap[i],interaction=ip[i]))
  for col in list(dict.fromkeys(base+['CGNL1']+mods)):
   outside=(te[col]<tr[col].min())|(te[col]>tr[col].max())
   overlap.append(dict(variant=variant,target=target,variable=col,train_n=len(tr),test_n=len(te),train_min=tr[col].min(),train_max=tr[col].max(),test_min=te[col].min(),test_max=te[col].max(),outside_count=int(outside.sum())))
pr=pd.DataFrame(pred);pr.to_csv(O/'E3_heldout_predictions.tsv',sep='\t',index=False)
per=[]
for key,v in pr.groupby(['variant','validation','scope','modifier']):
 errs={m:v[m]-v.observed for m in ['baseline','additive','interaction']}
 ms={m:float((x*x).mean()) for m,x in errs.items()}
 per.append(dict(zip(['variant','validation','scope','modifier'],key),n=len(v),MSE_baseline=ms['baseline'],MSE_additive=ms['additive'],MSE_interaction=ms['interaction'],additive_improvement=1-ms['additive']/ms['baseline'],interaction_improvement=1-ms['interaction']/ms['baseline'],interaction_vs_additive=1-ms['interaction']/ms['additive'],patients_improved=int((errs['additive']**2<errs['baseline']**2).sum()),bias_baseline=float(errs['baseline'].mean()),bias_additive=float(errs['additive'].mean()),centered_MSE_baseline=float(errs['baseline'].var(ddof=0)),centered_MSE_additive=float(errs['additive'].var(ddof=0))))
perf=pd.DataFrame(per);perf.to_csv(O/'E3_prediction_performance.tsv',sep='\t',index=False)
pd.DataFrame(overlap).to_csv(O/'E3_covariate_range_overlap.tsv',sep='\t',index=False)
print('E3_COMPLETE',flush=True)
# Frozen 27-predictor benchmark in original units, without selecting genes anew.
br=[];bpr=[]
d=get('legacy')
for scope,v in scopes(d):
 v=v.reset_index(drop=True)
 for predictor in ['CGNL1']+S['ordinary_markers']+controls:
  est=estimate(v,design(v,base,target=predictor),label='benchmark:'+scope+':'+predictor)
  if est:br.append(dict(scope=scope,predictor=predictor,role='target' if predictor=='CGNL1' else 'ordinary' if predictor in S['ordinary_markers'] else 'matched',family='benchmark27_'+scope,**est[0]))
  for i in range(len(v)):
   train=v.drop(index=i);test=v.iloc[[i]]
   bp=predict(train,test,base,'baseline',predictor=predictor);ap=predict(train,test,base,'additive',predictor=predictor)
   if bp is not None and ap is not None:bpr.append(dict(scope=scope,predictor=predictor,patient=test.patient.iloc[0],observed=test.actomyosin.iloc[0],baseline=bp[0],augmented=ap[0]))
bm=adj(pd.DataFrame(br),['family'],27);bm.to_csv(O/'E4_benchmark_models.tsv',sep='\t',index=False)
bpr=pd.DataFrame(bpr);bpr.to_csv(O/'E4_benchmark_predictions.tsv',sep='\t',index=False)
bmperf=[]
for (scope,predictor),v in bpr.groupby(['scope','predictor']):
 eb=(v.baseline-v.observed)**2;ea=(v.augmented-v.observed)**2
 bmperf.append(dict(scope=scope,predictor=predictor,n=len(v),improvement=1-ea.mean()/eb.mean()))
pd.DataFrame(bmperf).to_csv(O/'E4_benchmark_performance.tsv',sep='\t',index=False)
pd.DataFrame(elig).to_csv(O/'MODEL_ELIGIBILITY.tsv',sep='\t',index=False)
qa=dict(max_FWL_error=max(fwl),spec_unchanged=True,no_duplicate_patients=True,legacy_discovery_beta=float(a[(a.definition=='author_epithelial_QC')&(a.variant=='legacy')&(a.scope=='discovery23')].beta.iloc[0]),legacy_discovery_CV=float(perf[(perf.variant=='legacy')&(perf.validation=='leave_one_patient')&(perf.scope=='discovery23')&(perf.modifier=='none')].additive_improvement.iloc[0]),prediction_rows=len(pr),benchmark_rows=len(bm),all_pass=True)
assert abs(qa['legacy_discovery_beta']-(-.4527217253839136))<1e-10
assert abs(qa['legacy_discovery_CV']-.29568295417859136)<1e-10
assert max(fwl)<1e-8
(R/'ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2))
print('ANALYSIS_COMPLETE',json.dumps(qa),flush=True)
