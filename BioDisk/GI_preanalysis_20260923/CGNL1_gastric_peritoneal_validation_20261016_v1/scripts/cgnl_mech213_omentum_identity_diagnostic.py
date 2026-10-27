
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,itertools
import pandas as pd,numpy as np
from scipy import stats
p=Path.cwd();m=p/'mechanism_investigation_v2';out=m/'results/GSE286057_comparator_198'
qa=json.loads((out/'RAW_MATRIX_INDEPENDENT_QA_211.json').read_text(encoding='utf-8'));assert qa['pass'] and qa['samples']==36
spec=json.loads((m/'GSE286057_IDENTITY_DIAGNOSTIC_SPECIFICATION_213.json').read_text(encoding='utf-8'))
src=pd.read_csv(out/'source_profiles.tsv',sep=chr(9));ctx=pd.read_csv(out/'context_genes.tsv',sep=chr(9))
assert not ctx[['title','scheme','group','gene']].duplicated().any()
wide=ctx.pivot(index=['title','scheme','group'],columns='gene',values='log2CPM1').reset_index()
data=src.merge(wide,on=['title','scheme','group'],how='inner',validate='one_to_one')
data=data[data.n_cells>=20].copy();assert np.allclose(data.log2CPM1,data.CGNL1)
for name,genes in spec['panels'].items():
 assert data[genes].notna().all().all()
 data[name+'_comparator']=data.CGNL1-data[genes].mean(axis=1)
rows=[];samplerows=[]
for scheme in ['primary','sensitivity']:
 for group in ['strict_mesothelial','mesothelial_both_specific','mesothelial_twoUMI','mesothelial_no_epithelial_flag']:
  for population in ['all_deposited','exclude_borderline_and_relapse_cases','primary_high_grade_vs_benign_only']:
   d=data[(data.scheme==scheme)&(data.group==group)].copy()
   a=d[(d.site=='peritumoral')&(d.cohort=='ovarian_metastatic')]
   b=d[(d.site=='proximal')&(d.cohort=='nonmetastatic')]
   if population!='all_deposited':a=a[~a.donor.isin(['P9','P10'])]
   if population=='primary_high_grade_vs_benign_only':b=b[b.donor.isin(['P1','P2','P3','P4'])]
   assert not a.donor.duplicated().any() and not b.donor.duplicated().any()
   for measure in ['CGNL1','five_identity_comparator','two_specific_comparator']:
    av=a[measure].dropna().to_numpy();bv=b[measure].dropna().to_numpy()
    r={'scheme':scheme,'group':group,'population':population,'measure':measure,'n_cases':len(av),'n_controls':len(bv),'mean_case':float(av.mean()) if len(av) else None,'mean_control':float(bv.mean()) if len(bv) else None,'mean_difference':float(av.mean()-bv.mean()) if len(av) and len(bv) else None,'P':None,'FDR_two_comparators':None,'scope':'posthoc_descriptive'}
    if scheme=='primary' and group=='strict_mesothelial' and population=='all_deposited' and measure!='CGNL1' and min(len(av),len(bv))>=3:
     test=stats.permutation_test((av,bv),lambda a,b:np.mean(a)-np.mean(b),n_resamples=np.inf,alternative='two-sided',permutation_type='independent',vectorized=False)
     r.update(P=float(test.pvalue),scope='posthoc_exploratory_exact_permutation_family2')
    rows.append(r)
res=pd.DataFrame(rows);mask=res.P.notna();assert mask.sum()==2
res.loc[mask,'FDR_two_comparators']=stats.false_discovery_control(res.loc[mask,'P'].to_numpy())
data.to_csv(out/'identity_sample_values_213.tsv',sep=chr(9),index=False,encoding='utf-8')
res.to_csv(out/'identity_and_clinical_sensitivity_213.tsv',sep=chr(9),index=False,encoding='utf-8')
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'diagnostic_complete_not_final_mainline','raw_matrix_QA_passed':True,'primary_full_cohort':res[(res.scheme=='primary')&(res.group=='strict_mesothelial')&(res.population=='all_deposited')].replace({np.nan:None}).to_dict('records'),'primary_clinical_sensitivity':res[(res.scheme=='primary')&(res.group=='strict_mesothelial')&(res.population!='all_deposited')&(res.measure=='CGNL1')].replace({np.nan:None}).to_dict('records'),'causal_mechanism_established':False}
(out/'IDENTITY_DIAGNOSTIC_COMPLETE_213.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
