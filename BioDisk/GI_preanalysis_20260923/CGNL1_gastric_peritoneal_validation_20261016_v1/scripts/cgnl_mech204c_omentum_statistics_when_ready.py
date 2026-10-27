
from pathlib import Path
import subprocess
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,select,itertools,math
import numpy as np,pandas as pd
from scipy import stats
p=Path.cwd();m=p/'mechanism_investigation_v2'
launch=json.loads((m/'metadata/GSE286057_ANALYSIS_LAUNCH_198c.json').read_text(encoding='utf-8'));receipt=Path(launch['receipt'])
if not receipt.exists():
 try:
  subprocess.run(['/usr/bin/python3',str(p/'scripts/cgnl_wait_pidfd_system.py'),str(launch['wrapper_pid'])],check=True)
 except ProcessLookupError:pass
assert receipt.exists() and json.loads(receipt.read_text(encoding='utf-8'))['exit_code']==0
out=m/'results/GSE286057_comparator_198';assert json.loads((out/'SOURCE_COMPLETE.json').read_text())['all_pass']
x=pd.read_csv(out/'source_profiles.tsv',sep='	');qc=pd.read_csv(out/'QC.tsv',sep='	')
assert x[['title','scheme','group']].duplicated().sum()==0
# Independent validation of outputs against saved per-cell measurements.
checks=[]
for f in sorted((m/'processed/GSE286057_comparator_198').glob('*_cells.tsv.gz')):
 c=pd.read_csv(f,sep='	');prefix=f.name[:-len('_cells.tsv.gz')];acc=prefix.split('_',1)[0];s=x[x.accession==acc]
 assert len(s)==22 and not c.barcode.duplicated().any()
 maxerr=0.
 for _,r in s.iterrows():
  keep=c['qc_'+r.scheme].astype(bool);mask=c[r.group].astype(bool)&keep;n=int(mask.sum());tu=int(c.loc[mask,'CGNL1'].sum());total=int(c.loc[mask,'UMI'].sum())
  assert n==r.n_cells and tu==r.target_UMI and total==r.total_UMI
  if total:err=abs(math.log2(1+tu/total*1e6)-r.log2CPM1);assert err<1e-10;maxerr=max(maxerr,err)
  assert int(c.loc[mask,'epithelial_flag'].sum())==r.epithelial_flag_cells
 checks.append({'sample':s.title.iloc[0],'source_rows_checked':len(s),'max_log2CPM1_error':maxerr})
assert len(checks)==36
summary=[];pairs=[]
for scheme in ['primary','sensitivity']:
 for group in ['strict_mesothelial','mesothelial_both_specific','mesothelial_twoUMI','mesothelial_no_epithelial_flag']:
  for population in ['all_deposited_metastatic','primary_high_grade_only']:
   d=x[(x.scheme==scheme)&(x.group==group)&(x.n_cells>=20)&x.site.isin(['peritumoral','distant'])].copy()
   if population=='primary_high_grade_only':d=d[~d.donor.isin(['P9','P10'])]
   w=d.pivot(index='donor',columns='site',values='log2CPM1')
   good=w.dropna(subset=['peritumoral','distant']) if set(['peritumoral','distant']).issubset(w.columns) else pd.DataFrame()
   delta=(good.peritumoral-good.distant) if len(good) else pd.Series(dtype=float)
   primary=scheme=='primary' and group=='strict_mesothelial' and population=='all_deposited_metastatic'
   rec={'scheme':scheme,'group':group,'population':population,'role':'prespecified_primary' if primary else 'descriptive_sensitivity','n_paired_donors':len(delta),'mean_delta_log2CPM1':float(delta.mean()) if len(delta) else None,'median_delta_log2CPM1':float(delta.median()) if len(delta) else None,'positive':int((delta>0).sum()),'negative':int((delta<0).sum()),'P':None,'minimum_attainable_P':2**(1-len(delta)) if len(delta) else None,'epithelial_flags_retained':True}
   if primary and len(delta)>=3 and (delta!=0).any():
    nz=delta[delta!=0].to_numpy();ranks=stats.rankdata(abs(nz));obs=abs(np.dot(np.sign(nz),ranks))
    rec['P']=float(np.mean([abs(np.dot(sign,ranks))>=obs-1e-12 for sign in itertools.product([-1,1],repeat=len(nz))]))
   for donor,value in delta.items():
    rr=d[d.donor==donor].set_index('site')
    pairs.append({'scheme':scheme,'group':group,'population':population,'donor':donor,'delta_log2CPM1':float(value),'peritumoral_cells':int(rr.loc['peritumoral','n_cells']),'distant_cells':int(rr.loc['distant','n_cells']),'peritumoral_epithelial_flag_fraction':float(rr.loc['peritumoral','epithelial_flag_fraction']),'distant_epithelial_flag_fraction':float(rr.loc['distant','epithelial_flag_fraction'])})
   summary.append(rec)
# Primary secondary comparisons remain separately defined families.
base=x[(x.scheme=='primary')&(x.group=='strict_mesothelial')];secondary=[]
for measure,family in [('log2CPM1','within_source_expression'),('cell_fraction','mesothelial_proportion')]:
 d=base[base.n_cells>=20] if measure=='log2CPM1' else base
 a=d[(d.site=='peritumoral')&(d.cohort=='ovarian_metastatic')];b=d[(d.site=='proximal')&(d.cohort=='nonmetastatic')]
 assert not a.donor.duplicated().any() and not b.donor.duplicated().any()
 av=a[measure].dropna().to_numpy();bv=b[measure].dropna().to_numpy()
 rec={'measure':measure,'family':family,'n_metastatic_donors':len(av),'n_nonmetastatic_donors':len(bv),'mean_difference':float(av.mean()-bv.mean()) if len(av) and len(bv) else None,'P':None,'FDR':None,'test':'two-sided exact permutation of donor-level mean difference','scope':'heterogeneous diagnoses and region sampling; not all healthy controls; marker-defined'}
 if min(len(av),len(bv))>=3 and np.ptp(np.r_[av,bv])>0:
  res=stats.permutation_test((av,bv),lambda a,b:np.mean(a)-np.mean(b),permutation_type='independent',alternative='two-sided',n_resamples=np.inf,vectorized=False)
  rec['P']=float(res.pvalue);rec['FDR']=float(res.pvalue)
 secondary.append(rec)
pd.DataFrame(summary).to_csv(out/'paired_summary_204.tsv',sep='	',index=False,encoding='utf-8')
pd.DataFrame(pairs).to_csv(out/'paired_donors_204.tsv',sep='	',index=False,encoding='utf-8')
pd.DataFrame(secondary).to_csv(out/'secondary_contrasts_204.tsv',sep='	',index=False,encoding='utf-8')
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'scope':'all36 saved cell tables vs22 source rows per sample; original matrix integrity verified by downloader/source reader, independent raw-triplet validation still pending','checks':checks}
(out/'AGGREGATE_QA_204.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
(out/'STATISTICS_COMPLETE_204.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary':summary[0],'secondary':secondary,'original_matrix_independent_QA_pending':True,'causal_validation':False},indent=2),encoding='utf-8')
print(json.dumps({'paired':summary,'secondary':secondary},indent=2),flush=True)
