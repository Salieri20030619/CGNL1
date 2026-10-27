from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,itertools,datetime
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');O=R/'results/platform_adapted';O.mkdir(exist_ok=True)
s=json.loads((R/'SPECIFICATION.json').read_text());panels=json.loads((R/'MEASUREMENT_AMENDMENT01.json').read_text())['common_measured_panels'];genes=['CGNL1']+s['ordinary_marker_comparators']
tables=[];regionrows=[]
for f in sorted((R/'cache').glob('C_*.npz')):
 z=np.load(f,allow_pickle=True);g=z['genes'].astype(str);umi=z['umi'];X=z['counts']
 sample=str(z['sample'][0]);patient=int(z['patient'][0])
 V=np.log1p(X/np.maximum(umi,1)[:,None]*1e4)
 d=pd.DataFrame({'sample':sample,'patient':patient,'region':z['region'],'logUMI':np.log1p(umi)})
 d['micro']=d.region.eq('Liver micrometastasis tumor').astype(float)
 for name,gs in panels.items():
  ix=[list(g).index(q) for q in gs if q in g]
  d[name]=V[:,ix].mean(1) if len(ix)>=.8*len(gs) else np.nan
 for q in genes:d[q]=V[:,list(g).index(q)] if q in g else np.nan
 tables.append(d)
 for region in np.unique(z['region']):
  ix=z['region']==region;cnt=X[ix].sum(0);den=umi[ix].sum();vals=np.log1p(cnt/max(den,1)*1e4)
  rr=dict(patient=patient,sample=sample,region=region,spots=int(ix.sum()),UMI=float(den))
  for q in genes:
   if q in g:
    j=list(g).index(q);rr[q]=float(vals[j]);rr[q+'_count']=float(cnt[j]);rr[q+'_detection']=float((X[ix,j]>0).mean())
  for name,gs in panels.items():
   gi=[list(g).index(q) for q in gs if q in g];rr[name]=float(vals[gi].mean()) if len(gi)>=.8*len(gs) else np.nan
  regionrows.append(rr)
d=pd.concat(tables,ignore_index=True);pd.DataFrame(regionrows).to_csv(O/'C_sample_region_profiles.tsv',sep='\t',index=False)
d.to_csv(O/'C_spot_model_inputs.tsv.gz',sep='\t',index=False)
rows=[];eligible=[];numeric=[]
for threshold in [10,20]:
 for patient,v in d.groupby('patient'):
  nm=int(v.micro.sum());na=len(v)-nm
  for model,covs in s['C']['models'].items():
   info=dict(patient=patient,model=model,min_micro=threshold,micro_spots=nm,macro_spots=na)
   if nm<threshold or na<20:
    eligible.append(dict(info,status='insufficient_pathology_region_coverage'));continue
   v=v.reset_index(drop=True);Z=pd.DataFrame(index=v.index)
   for c in covs:
    if c=='section':
     for k in sorted(v['sample'].unique())[1:]:Z['section_'+k]=v['sample'].eq(k).astype(float)
    else:Z[c]=v[c]
   if not np.isfinite(Z.to_numpy()).all():
    eligible.append(dict(info,status='insufficient_panel_coverage'));continue
   raw=Z.to_numpy(dtype=float);sd=raw.std(0,ddof=1);active=sd>1e-10
   F=np.c_[np.ones(len(v)),(raw[:,active]-raw[:,active].mean(0))/sd[active],v.micro.to_numpy()]
   rank=np.linalg.matrix_rank(F);cond=np.linalg.cond(F)
   if rank<F.shape[1] or len(v)-F.shape[1]<10 or cond>1e8:
    eligible.append(dict(info,status='stage_not_identifiable_after_covariates',rank=rank,parameters=F.shape[1],condition_number=float(cond)));continue
   eligible.append(dict(info,status='estimated',rank=rank,parameters=F.shape[1],condition_number=float(cond)))
   for gene in genes:
    y=v[gene].to_numpy();ys=y.std(ddof=1)
    if not np.isfinite(y).all() or ys<1e-10:continue
    fit=np.linalg.lstsq(F,y,rcond=None)[0];co=float(fit[-1])
    rows.append(dict(info,gene=gene,beta_logCP10k=co,beta_within_patient_SD=co/ys,target_nonzero_spots=int((y>0).sum()),condition_number=float(cond)))
    if gene=='CGNL1':
     N=F[:,:-1];rx=F[:,-1]-N@np.linalg.lstsq(N,F[:,-1],rcond=None)[0];ry=y-N@np.linalg.lstsq(N,y,rcond=None)[0]
     numeric.append(dict(info,FWL_error=abs(float(rx@ry/(rx@rx))-co)))
q=pd.DataFrame(rows);q.to_csv(O/'C_patient_coefficients.tsv',sep='\t',index=False)
pd.DataFrame(eligible).to_csv(O/'C_model_eligibility.tsv',sep='\t',index=False)
def planned_bh3(values):
 values=np.asarray(values,dtype=float)
 return multipletests(np.r_[values,np.ones(max(0,3-len(values)))],method='fdr_bh')[1][:len(values)]
def signflip(x):
 x=np.asarray(x,dtype=float);n=len(x)
 if n<2:return np.nan
 observed=abs(x.mean());signs=np.array(list(itertools.product([-1,1],repeat=n)))
 return float((abs(signs@x/n)>=observed-1e-12).mean())
def summarise(x):
 x=np.asarray(x,dtype=float);n=len(x)
 loo=[np.delete(x,i).mean() for i in range(n)] if n>1 else [np.nan]
 return dict(n=n,mean_effect=float(x.mean()),median_effect=float(np.median(x)),positive=int((x>0).sum()),negative=int((x<0).sum()),P_signflip=signflip(x),P_ttest=float(stats.ttest_1samp(x,0).pvalue) if n>1 and np.std(x)>0 else np.nan,LOO_mean_min=float(np.nanmin(loo)),LOO_mean_max=float(np.nanmax(loo)))
summ=[]
for (threshold,model,gene),v in q.groupby(['min_micro','model','gene']):
 summ.append(dict(min_micro=threshold,model=model,gene=gene,scale='ln1pCP10k',**summarise(v.beta_logCP10k)))
t=pd.DataFrame(summ)
for method in ['signflip','ttest']:
 t['FDR_benchmark_'+method]=np.nan;t['FDR_target_'+method]=np.nan
 for key,ix in t.groupby(['min_micro','model']).groups.items():
  good=[i for i in ix if np.isfinite(t.loc[i,'P_'+method])]
  if good:t.loc[good,'FDR_benchmark_'+method]=multipletests(t.loc[good,'P_'+method],method='fdr_bh')[1]
 for threshold,ix in t[t.gene.eq('CGNL1')].groupby('min_micro').groups.items():
  good=[i for i in ix if np.isfinite(t.loc[i,'P_'+method])]
  if good:t.loc[good,'FDR_target_'+method]=planned_bh3(t.loc[good,'P_'+method])
t.to_csv(O/'C_gene_summary.tsv',sep='\t',index=False)
special=[];spatients=[]
for (threshold,model),v in q.groupby(['min_micro','model']):
 mat=v.pivot(index='patient',columns='gene',values='beta_within_patient_SD').reindex(columns=genes).dropna()
 if mat.empty:continue
 diff=mat.CGNL1-mat[s['ordinary_marker_comparators']].median(axis=1)
 for patient,value in diff.items():spatients.append(dict(min_micro=threshold,model=model,patient=patient,standardized_target_minus_median_marker=value))
 special.append(dict(min_micro=threshold,model=model,scale='within_patient_SD_target_minus_median_marker',**summarise(diff)))
st=pd.DataFrame(special)
for method in ['signflip','ttest']:
 st['FDR_'+method]=np.nan
 for _,ix in st.groupby('min_micro').groups.items():
  good=[i for i in ix if np.isfinite(st.loc[i,'P_'+method])]
  if good:st.loc[good,'FDR_'+method]=planned_bh3(st.loc[good,'P_'+method])
st.to_csv(O/'C_specificity_summary.tsv',sep='\t',index=False);pd.DataFrame(spatients).to_csv(O/'C_specificity_patients.tsv',sep='\t',index=False)
pd.DataFrame(numeric).to_csv(O/'C_FWL_numeric_check.tsv',sep='\t',index=False)
assert max(x['FWL_error'] for x in numeric)<1e-8
(O/'C_ANALYSIS_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'patient_coefficient_rows':len(q),'FWL_max_error':max(x['FWL_error'] for x in numeric),'biological_unit':'patient','spot_inferential_P_used':False},indent=2))
print('C_ANALYSIS_COMPLETE',flush=True)
