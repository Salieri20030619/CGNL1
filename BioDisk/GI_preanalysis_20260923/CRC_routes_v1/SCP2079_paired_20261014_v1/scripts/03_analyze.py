from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib,datetime
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent.parent;O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text(encoding='utf-8'))
assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text(encoding='utf-8'))['sha256']
D=pd.read_csv(O/'patient_deltas.tsv',sep='\t')
models=[];predictions=[];performance=[];changes=[];descriptive=[]
def populations(df):
 return {'fixed_qc_ge20_all_sites':df[df.definition.eq('fixed_qc')&df.min_cells.ge(20)],'fixed_qc_ge50_all_sites':df[df.definition.eq('fixed_qc')&df.min_cells.ge(50)],'author_retained_ge20_all_sites':df[df.definition.eq('author_retained')&df.min_cells.ge(20)],'fixed_qc_ge20_liver':df[df.definition.eq('fixed_qc')&df.min_cells.ge(20)&df.tissueSite.eq('Liver')]}
def design(train,test,features):
 mu=train[features].mean();sd=train[features].std(ddof=1)
 if (sd<1e-12).any() or not np.isfinite(sd).all():return None
 X=np.column_stack([np.ones(len(train)),(train[features]-mu)/sd]);T=np.column_stack([np.ones(len(test)),(test[features]-mu)/sd])
 if not np.isfinite(X).all() or not np.isfinite(T).all() or np.linalg.matrix_rank(X)<X.shape[1] or np.linalg.cond(X)>1e8:return None
 return X,T,sd

def signed_rank_exact(delta):
 d=np.asarray(delta);d=d[d!=0]
 if len(d)==0:return 1.,0,0
 ranks=np.rint(stats.rankdata(np.abs(d))*2).astype(int)
 # Exact conditional distribution via subset-sum dynamic programming.
 counts=np.array([1],dtype=np.int64)
 for w in ranks:
  new=np.zeros(len(counts)+w,dtype=np.int64);new[:len(counts)]+=counts;new[w:]+=counts;counts=new
 observed=int(ranks[d>0].sum());p=min(1.,2*min(counts[:observed+1].sum(),counts[observed:].sum())/2**len(d))
 return float(p),float(min(observed,ranks.sum()-observed)/2),len(d)
for measurement,df in D.groupby('measurement'):
 for population,g in populations(df).items():
  g=g.reset_index(drop=True);n=len(g)
  assert g.donor_id.is_unique
  for outcome in S['outcomes']:
   rho=float(stats.spearmanr(g.CGNL1,g[outcome]).statistic) if n>1 and g.CGNL1.std()>0 and g[outcome].std()>0 else np.nan
   descriptive.append({'measurement':measurement,'population':population,'outcome':outcome,'n':n,'spearman_rho_descriptive':rho,'CGNL1_positive_changes':int((g.CGNL1>0).sum()),'outcome_positive_changes':int((g[outcome]>0).sum())})
   for model,controls in S['models'].items():
    features=controls+['CGNL1'];mat=design(g,g,features) if n>1 else None
    family=measurement+'|'+population+'|'+model+'|CGNL1_two_outcomes'
    info={'measurement':measurement,'population':population,'model':model,'outcome':outcome,'n':n,'family':family}
    row={**info,'status':'not_estimable','beta_standardized':np.nan,'slope':np.nan,'CI_low':np.nan,'CI_high':np.nan,'P':np.nan,'FDR':np.nan,'df':np.nan}
    eligible=n>=S['inferential_min_n'] and mat is not None and n-mat[0].shape[1]>=S['inferential_min_residual_df'] and g[outcome].std()>1e-12
    if not eligible:
     row['status']='descriptive_only_n_lt_10' if n<10 else 'rank_df_or_variance_failure';models.append(row);continue
    X,T,sd=mat;fit=sm.OLS(g[outcome].to_numpy(),X).fit(cov_type='HC3',use_t=True);i=len(features);ci=fit.conf_int()[i]/sd.CGNL1
    row.update(status='estimated',beta_standardized=float(fit.params[i]/g[outcome].std()),slope=float(fit.params[i]/sd.CGNL1),CI_low=float(ci[0]),CI_high=float(ci[1]),P=float(fit.pvalues[i]),df=float(fit.df_resid));models.append(row)
    fold=[]
    for held in range(n):
     tr=g.drop(index=held);te=g.iloc[[held]];base=design(tr,te,controls);aug=design(tr,te,features)
     if base is None or aug is None or len(tr)-aug[0].shape[1]<S['prediction_min_training_residual_df']:continue
     bp=float((base[1]@np.linalg.lstsq(base[0],tr[outcome],rcond=None)[0])[0]);ap=float((aug[1]@np.linalg.lstsq(aug[0],tr[outcome],rcond=None)[0])[0])
     fold.append({**info,'donor_id':te.donor_id.iloc[0],'observed':float(te[outcome].iloc[0]),'baseline':bp,'augmented':ap})
    predictions.extend(fold)
    perf={**info,'completed_folds':len(fold),'all_folds_complete':len(fold)==n,'MSE_baseline':np.nan,'MSE_augmented':np.nan,'relative_MSE_improvement':np.nan,'patients_improved':np.nan}
    if len(fold)==n:
     f=pd.DataFrame(fold);be=(f.baseline-f.observed)**2;ae=(f.augmented-f.observed)**2
     perf.update(MSE_baseline=float(be.mean()),MSE_augmented=float(ae.mean()),relative_MSE_improvement=float(1-ae.mean()/be.mean()),patients_improved=int((ae<be).sum()))
    performance.append(perf)
  for readout in ['CGNL1']+list(S['outcomes']):
   delta=g[readout];p,stat,nz=signed_rank_exact(delta) if n>=10 else (np.nan,np.nan,int((delta!=0).sum()))
   changes.append({'measurement':measurement,'population':population,'readout':readout,'n':n,'mean_delta':float(delta.mean()),'median_delta':float(delta.median()),'positive_patients':int((delta>0).sum()),'negative_patients':int((delta<0).sum()),'zero_patients':int((delta==0).sum()),'P':p,'FDR':np.nan,'statistic':stat,'nonzero_pairs':nz,'status':'tested' if n>=10 else 'descriptive_only_n_lt_10','family':measurement+'|'+population+'|three_paired_changes'})
M=pd.DataFrame(models);C=pd.DataFrame(changes)
for table,size in [(M,2),(C,3)]:
 for family,ids in table.groupby('family').groups.items():
  valid=table.loc[ids,'P'].dropna();p=np.r_[valid.to_numpy(),np.ones(size-len(valid))]
  if len(valid):table.loc[valid.index,'FDR']=multipletests(p,method='fdr_bh')[1][:len(valid)]
 table['nominal_P_lt_0_05']=table.P<.05;table['FDR_lt_0_05']=table.FDR<.05
M.to_csv(O/'CGNL1_delta_models.tsv',sep='\t',index=False,encoding='utf-8')
C.to_csv(O/'paired_changes.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(predictions).to_csv(O/'heldout_predictions.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(performance).to_csv(O/'predictive_increment.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(descriptive).to_csv(O/'descriptive_associations.tsv',sep='\t',index=False,encoding='utf-8')
primary_measure='pooled_raw_CP10k' if S['scale']=='verified_raw_counts' else 'mean_author_log'
primary=M[M.measurement.eq(primary_measure)&M.population.eq('fixed_qc_ge20_all_sites')&M.model.eq('M1')]
print(primary.to_string(index=False),flush=True)
print(C[C.measurement.eq(primary_measure)&C.population.eq('fixed_qc_ge20_all_sites')].to_string(index=False),flush=True)
complete={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':len(M),'estimated_models':int(M.status.eq('estimated').sum()),'prediction_rows':len(predictions),'all_predictive_models_have_complete_folds':all(x['all_folds_complete'] for x in performance),'primary_models':json.loads(primary.to_json(orient='records'))}
(R/'ANALYSIS_COMPLETE.json').write_text(json.dumps(complete,indent=2),encoding='utf-8')