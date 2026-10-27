from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime
import pandas as pd,numpy as np
from scipy import stats
import statsmodels.api as sm
O=Path.cwd()/'mechanism_investigation_v2';B=O.parent.parent;panels=json.loads((O.parent/'IMPLEMENTATION_FREEZE.json').read_text())['bulk_panels'];spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'issue':'GC-Q06','purpose':'sensitivitytocaseGC05samplingchoice,notnewindependentcohort','assumption':'sharedGC05patientgrouping inheritedfromoldmetadata; clinicaltimepointconfirmationpending','variants':['original_a','replace_with_b','average_a_b_within_patient','exclude_GC05'],'retain_independent_n':'17or16,nottreat18samplepairsas18patients','tests':'pairedexactrankandpreservationadjustedHC3targetintercept; BHacross4variantsseparatelybytestfamily','interpretation':'robustnessauditonly; do notselectbestP orreplaceoriginalfamilies/results'};(O/'GC05_SENSITIVITY_SPECIFICATION.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
x=pd.read_csv(B/'03_processed/GC314812_TMM_logCPM.tsv.gz',sep='\t');q=x[x.symbol=='CGNL1'];assert len(q)==1;meta=pd.read_csv(B/'01_metadata/GSE314812_sample_design.tsv',sep='\t');d=pd.read_csv(O.parent/'results/bulk_patient_deltas.tsv',sep='\t');assert len(d)==17 and (d.patient=='GC05').sum()==1
changes={};values=[]
for suffix in ['a','b']:
 pm='GC05'+suffix+'_GCPM';pt='GC05'+suffix+'_GCPT';r={'CGNL1':float(q[pm].iloc[0]-q[pt].iloc[0])}
 for panel,gs in panels.items():
  z=x[x.symbol.isin(gs)].groupby('symbol')[[pm,pt]].mean();r[panel]=float(z[pm].mean()-z[pt].mean())
 pp=meta.set_index('sample');r['preservation_delta']=int(pp.loc[pm,'preservation method']=='FFPE')-int(pp.loc[pt,'preservation method']=='FFPE');changes[suffix]=r;values.append(dict(sampling=suffix,primary=float(q[pt].iloc[0]),metastasis=float(q[pm].iloc[0]),**r))
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
r.to_csv(O/'results/GSE314812_GC05_sensitivity.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(values).to_csv(O/'results/GSE314812_GC05_sampling_values.tsv',sep='\t',index=False,encoding='utf-8');print(r.to_string(index=False));print('GC05',pd.DataFrame(values).to_string(index=False))
