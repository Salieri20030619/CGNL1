from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,itertools,datetime,hashlib
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');O=R/'results'
spec=json.loads((R/'SPECIFICATION.json').read_text())
def adjust(t,group,pcol,outcol):
 t[outcol]=np.nan
 for _,ix in t.groupby(group,dropna=False).groups.items():
  ok=[i for i in ix if np.isfinite(t.loc[i,pcol])]
  if ok:t.loc[ok,outcol]=multipletests(t.loc[ok,pcol],method='fdr_bh')[1]
 return t
def design(v,cols,predictor=None):
 d=pd.DataFrame(index=v.index)
 for c in cols:
  if c=='study':
   for name in sorted(v.dataset.unique())[1:]:d['study_'+name]=v.dataset.eq(name).astype(float)
  else:d[c]=v[c].astype(float)
 if predictor:d['added']=v[predictor].astype(float)
 return d
def fit_arrays(X,y,Xtest=None):
 # All preprocessing depends only on training X.
 X=np.asarray(X,dtype=float);y=np.asarray(y,dtype=float)
 mu=X.mean(0);sd=X.std(0,ddof=1);active=np.isfinite(sd)&(sd>1e-10)
 z=(X[:,active]-mu[active])/sd[active];F=np.c_[np.ones(len(y)),z]
 if np.linalg.matrix_rank(F)<F.shape[1] or len(y)-F.shape[1]<3 or np.linalg.cond(F)>1e8:return None
 fit=np.linalg.lstsq(F,y,rcond=None)[0]
 if Xtest is None:return F,fit,active,mu,sd
 zt=(np.asarray(Xtest,dtype=float)[:,active]-mu[active])/sd[active]
 return np.c_[np.ones(len(zt)),zt]@fit
pt=pd.read_csv(O/'A_patient_profiles.tsv',sep='\t')
gene=pd.read_csv(O/'A_patient_allgenes_logCP10k.tsv.gz',sep='\t',index_col=0)
controls=json.loads((O/'A_CONTROL_FREEZE.json').read_text())['controls']
predictors=['CGNL1']+spec['ordinary_marker_comparators']+controls
for p in predictors:
 if p not in pt:pt[p]=pt.patient.map(gene[p]) if p in gene else np.nan
pt.to_csv(O/'A_model_input.tsv',sep='\t',index=False)
rows=[];elig=[];pred=[];transfer=[];descriptive=[];numeric=[]
for scope in spec['A']['scopes']:
 v=pt.copy() if scope=='pooled' else pt[pt.dataset.eq(scope)].copy()
 v=v.reset_index(drop=True)
 for mn,cols in spec['A']['models'].items():
  for out in spec['outcomes']:
   base=design(v,cols)
   for p in predictors:
    role='CGNL1' if p=='CGNL1' else 'ordinary_marker' if p in spec['ordinary_marker_comparators'] else 'matched_expression'
    vv=v.dropna(subset=[out,p]+[c for c in cols if c!='study']).reset_index(drop=True)
    raw=design(vv,cols,p);y=vv[out].to_numpy();n=len(vv)
    status=dict(scope=scope,model=mn,outcome=out,predictor=p,n=n,role=role)
    res=fit_arrays(raw,y) if n>2 else None
    if res is None or not res[2][-1] or n-res[0].shape[1]<5 or np.std(y)<1e-10:
     elig.append(dict(status,status='not_estimable_rank_df_variation'));continue
    F,coef,active,mu,sd=res;df=n-F.shape[1];cond=float(np.linalg.cond(F))
    elig.append(dict(status,status='estimated',residual_df=df,condition_number=cond))
    for variance in ['HC3','OLS']:
     fit=sm.OLS(y,F).fit(cov_type='HC3' if variance=='HC3' else 'nonrobust',use_t=True)
     slope=float(fit.params[-1]/sd[-1]);ci=np.asarray(fit.conf_int())[-1]/sd[-1]
     rows.append(dict(status,variance=variance,slope=slope,beta_standardized=float(fit.params[-1]/np.std(y,ddof=1)),CI_low=float(ci[0]),CI_high=float(ci[1]),P=float(fit.pvalues[-1]),partial_R2=float(fit.tvalues[-1]**2/(fit.tvalues[-1]**2+df)) if variance=='OLS' else np.nan,residual_df=df,condition_number=cond))
    if p=='CGNL1':
     # Frisch-Waugh-Lovell verifies the target coefficient independently.
     N=F[:,:-1];rx=F[:,-1]-N@np.linalg.lstsq(N,F[:,-1],rcond=None)[0];ry=y-N@np.linalg.lstsq(N,y,rcond=None)[0]
     fwl=float(rx@ry/(rx@rx))
     numeric.append(dict(scope=scope,model=mn,outcome=out,abs_coefficient_error=abs(fwl-coef[-1])))
     nu=design(vv,cols).to_numpy();rr=fit_arrays(nu,vv.CGNL1.to_numpy())
     if rr:
      ff,bb,_,_,_=rr;yp=ff@bb;ss=float(((vv.CGNL1-vv.CGNL1.mean())**2).sum())
      descriptive.append(dict(scope=scope,model=mn,outcome=out,n=n,target_variance_explained=1-float(((vv.CGNL1-yp)**2).sum())/ss if ss>0 else np.nan))
    # Evaluate independent held-out patients; no held-out scaling or fitting.
    braw=design(vv,cols).to_numpy();araw=raw.to_numpy()
    for i in range(n):
     train=np.arange(n)!=i
     b=fit_arrays(braw[train],y[train],braw[[i]]);a=fit_arrays(araw[train],y[train],araw[[i]])
     if b is not None and a is not None:
      pred.append(dict(scope=scope,model=mn,outcome=out,predictor=p,role=role,patient=vv.patient.iloc[i],observed=float(y[i]),baseline=float(b[0]),augmented=float(a[0])))
print('A_MODELS_FITTED',len(rows),flush=True)
models=pd.DataFrame(rows)
models['family_target']='A_target_'+models.scope+'_'+models.model+'_'+models.variance
models['family_benchmark']='A_benchmark_'+models.scope+'_'+models.model+'_'+models.outcome+'_'+models.variance
models=adjust(models,['family_benchmark'],'P','FDR_benchmark')
models['FDR_target']=np.nan
for _,ix in models[models.predictor.eq('CGNL1')].groupby('family_target').groups.items():models.loc[ix,'FDR_target']=multipletests(models.loc[ix,'P'],method='fdr_bh')[1]
models.to_csv(O/'A_all_models.tsv',sep='\t',index=False);models[models.predictor.eq('CGNL1')].to_csv(O/'A_CGNL1_models.tsv',sep='\t',index=False)
pd.DataFrame(elig).to_csv(O/'A_model_eligibility.tsv',sep='\t',index=False)
pr=pd.DataFrame(pred);pr.to_csv(O/'A_heldout_predictions.tsv',sep='\t',index=False)
perf=[]
for key,g in pr.groupby(['scope','model','outcome','predictor','role']):
 b=(g.observed-g.baseline)**2;a=(g.observed-g.augmented)**2
 perf.append(dict(zip(['scope','model','outcome','predictor','role'],key),heldout_patients=len(g),MSE_baseline=float(b.mean()),MSE_augmented=float(a.mean()),relative_MSE_improvement=float(1-a.mean()/b.mean()) if b.mean()>0 else np.nan,patients_improved=int((a<b).sum()),median_loss_change=float((a-b).median())))
pd.DataFrame(perf).to_csv(O/'A_predictive_increment.tsv',sep='\t',index=False)
# Cohort transfer does not use target-cohort outcomes or target-cohort scaling.
cols=[x for x in spec['A']['models']['M1'] if x!='study']
for trainset,testset in [('GSE245552','GSE285991'),('GSE285991','GSE245552')]:
 tr=pt[pt.dataset.eq(trainset)];te=pt[pt.dataset.eq(testset)]
 for out in spec['outcomes']:
  for p in ['CGNL1']+spec['ordinary_marker_comparators']:
   tt=tr.dropna(subset=cols+[out,p]);ee=te.dropna(subset=cols+[out,p])
   if len(tt)<7 or len(ee)<2:continue
   b=fit_arrays(design(tt,cols),tt[out],design(ee,cols));a=fit_arrays(design(tt,cols,p),tt[out],design(ee,cols,p))
   if b is None or a is None:continue
   for i,(_,row) in enumerate(ee.iterrows()):transfer.append(dict(train_cohort=trainset,test_cohort=testset,outcome=out,predictor=p,patient=row.patient,observed=row[out],baseline=b[i],augmented=a[i]))
pd.DataFrame(transfer).to_csv(O/'A_cross_cohort_predictions.tsv',sep='\t',index=False)
pd.DataFrame(descriptive).drop_duplicates(['scope','model']).to_csv(O/'A_CGNL1_variance_explained.tsv',sep='\t',index=False)
pd.DataFrame(numeric).to_csv(O/'A_FWL_numeric_check.tsv',sep='\t',index=False)
assert max(x['abs_coefficient_error'] for x in numeric)<1e-8
# Paired disease-site profiles and fixed marker-relative contrasts.
bp=pd.read_csv(O/'B_all_group_profiles.tsv',sep='\t')
metrics=['CGNL1','CGNL1_minus_epithelial','CGNL1_minus_mature']
bp[metrics[1]]=bp.CGNL1-bp.epithelial_identity;bp[metrics[2]]=bp.CGNL1-bp.mature_colonocyte
pair=[]
for group in ['Tumor_no_Tu10','CNV95_noTu10','CNV99_noTu10']:
 v=bp[bp.celltype.eq(group)&bp.n_cells.ge(20)]
 for patient,g in v.groupby('patient'):
  if not {'CCT','LCT'}.issubset(set(g.site)):continue
  a=g[g.site.eq('CCT')].iloc[0];b=g[g.site.eq('LCT')].iloc[0]
  for metric in metrics:pair.append(dict(group=group,patient=patient,metric=metric,primary=a[metric],liver=b[metric],delta=b[metric]-a[metric],n_shared_states=np.nan))
for patient,v in bp[bp.celltype.str.match(r'^Tu[0-9]{2}_')&~bp.celltype.str.startswith('Tu10')&bp.n_cells.ge(20)].groupby('patient'):
 shared=set(v[v.site.eq('CCT')].celltype)&set(v[v.site.eq('LCT')].celltype)
 if len(shared)<2:continue
 for metric in metrics:
  a=v[v.site.eq('CCT')&v.celltype.isin(shared)][metric].mean();b=v[v.site.eq('LCT')&v.celltype.isin(shared)][metric].mean()
  pair.append(dict(group='shared_author_subtypes',patient=patient,metric=metric,primary=a,liver=b,delta=b-a,n_shared_states=len(shared)))
pa=pd.DataFrame(pair);pa.to_csv(O/'B_paired_effects.tsv',sep='\t',index=False)
summ=[]
for (group,metric),g in pa.groupby(['group','metric']):
 d=g.delta.to_numpy();d=d[np.isfinite(d)];nz=d[np.abs(d)>1e-12]
 summ.append(dict(group=group,metric=metric,n=len(d),positive=int((d>0).sum()),negative=int((d<0).sum()),mean_delta=float(np.mean(d)),median_delta=float(np.median(d)),P=stats.binomtest(int((nz>0).sum()),len(nz),.5).pvalue if len(nz)>=2 else np.nan,family='B_three_contrasts_'+group,evidence='small_sample_exploratory_site_contrast'))
summ=adjust(pd.DataFrame(summ),['family'],'P','FDR');summ.to_csv(O/'B_paired_summary.tsv',sep='\t',index=False)
(O/'AB_ANALYSIS_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'A_models':len(models),'A_heldout_predictions':len(pr),'B_paired_rows':len(pa),'FWL_max_error':max(x['abs_coefficient_error'] for x in numeric)},indent=2))
print('AB_ANALYSIS_COMPLETE',flush=True)
