from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime
import pandas as pd,numpy as np
from scipy import stats
import statsmodels.api as sm
O=Path.cwd()/'mechanism_investigation_v2';B=O.parent.parent;panels=json.loads((O.parent/'IMPLEMENTATION_FREEZE.json').read_text())['bulk_panels'];spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'issue':'GC-Q06','purpose':'sensitivitytocaseGC05samplingchoice,notnewindependentcohort','assumption':'sharedGC05patientgrouping inheritedfromoldmetadata; clinicaltimepointconfirmationpending','variants':['original_a','replace_with_b','average_a_b_within_patient','exclude_GC05'],'retain_independent_n':'17or16,nottreat18samplepairsas18patients','tests':'pairedexactrankandpreservationadjustedHC3targetintercept; BHacross4variantsseparatelybytestfamily','interpretation':'robustnessauditonly; do notselectbestP orreplaceoriginalfamilies/results'};assert (O/'GC05_SENSITIVITY_SPECIFICATION.json').exists()
amend={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'LegacyTMMmatrixexcludesGC05b. UsealreadydefinedrawcountCPMsensitivityforallvariants;do notmixTMMandCPM orreplaceprimaryresults.','measure':'log2(full-libraryCPM+0.5)','outcome_selection':False}
(O/'GC05_CPM_IMPLEMENTATION_AMENDMENT.json').write_text(json.dumps(amend,indent=2),encoding='utf-8')
meta=pd.read_csv(B/'01_metadata/GSE314812_sample_design.tsv',sep='\t');raw=pd.read_csv(B/'02_raw/GSE314812_counts.tsv.gz',sep='\t');q=raw[raw.gene_id.astype(str).str.split('.').str[0]=='ENSG00000128849'];assert len(q)==1;cols=meta['sample'].tolist();assert all(c in raw.columns for c in cols);cp=np.log2(q[cols].iloc[0]/raw[cols].sum()*1e6+.5)
d=pd.read_csv(O.parent/'results/bulk_patient_deltas.tsv',sep='\t');d['CGNL1']=[float(cp[r.PM_sample]-cp[r.PT_sample]) for _,r in d.iterrows()];assert len(d)==17 and (d.patient=='GC05').sum()==1;assert abs(d.CGNL1.mean()-json.loads((O.parent/'results/BULK_RESULT.json').read_text())['naive_CPM_mean_delta'])<1e-12
changes={};values=[]
for suffix in ['a','b']:
 pm='GC05'+suffix+'_GCPM';pt='GC05'+suffix+'_GCPT';r={'CGNL1':float(cp[pm]-cp[pt])}
 pp=meta.set_index('sample');r['preservation_delta']=int(pp.loc[pm,'preservation method']=='FFPE')-int(pp.loc[pt,'preservation method']=='FFPE');changes[suffix]=r;values.append(dict(sampling=suffix,primary=float(cp[pt]),metastasis=float(cp[pm]),**r))
assert abs(changes['a']['CGNL1']-d.loc[d.patient=='GC05','CGNL1'].iloc[0])<1e-12
rows=[]
for variant in spec['variants']:
 v=d.copy()
 if variant=='exclude_GC05':v=v[v.patient!='GC05'].copy()
 elif variant in ['replace_with_b','average_a_b_within_patient']:
  change=changes['b'] if variant=='replace_with_b' else {k:(changes['a'][k]+changes['b'][k])/2 for k in changes['a']}
  for k,val in change.items():v.loc[v.patient=='GC05',k]=val
 a=v.CGNL1.to_numpy();nz=a[a!=0];rank=stats.rankdata(abs(nz));dist=np.array([0.])
 for r in rank:dist=np.r_[dist,dist+r]
 obs=min(rank[nz>0].sum(),rank[nz<0].sum());P=min(1.,2*float(np.mean(dist<=obs+1e-12)));Z=np.c_[np.ones(len(v)),v.preservation_delta];fit=sm.OLS(a,Z).fit(cov_type='HC3',use_t=True)
 rows.append({'variant':variant,'n_patients_assumed':len(v),'mean_delta':a.mean(),'positive':int((a>0).sum()),'negative':int((a<0).sum()),'paired_rank_P':P,'preservation_adjusted_effect':float(fit.params[0]),'preservation_adjusted_P':float(fit.pvalues[0])})
r=pd.DataFrame(rows)
for p in ['paired_rank_P','preservation_adjusted_P']:r[p+'_FDR_four_variants']=stats.false_discovery_control(r[p])
r.to_csv(O/'results/GSE314812_GC05_CPM_sensitivity.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(values).to_csv(O/'results/GSE314812_GC05_CPM_sampling_values.tsv',sep='\t',index=False,encoding='utf-8');print(r.to_string(index=False));print('GC05',pd.DataFrame(values).to_string(index=False))
