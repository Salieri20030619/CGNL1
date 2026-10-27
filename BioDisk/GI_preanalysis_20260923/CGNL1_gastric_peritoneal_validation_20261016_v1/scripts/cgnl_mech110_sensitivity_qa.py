from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,itertools
import pandas as pd,numpy as np
from scipy.stats import wilcoxon,t
O=Path.cwd()/'mechanism_investigation_v2';B=O.parent.parent
raw=pd.read_csv(B/'02_raw/GSE314812_counts.tsv.gz',sep='\t');meta=pd.read_csv(B/'01_metadata/GSE314812_sample_design.tsv',sep='\t').set_index('sample');target=raw[raw.gene_id.str.split('.').str[0]=='ENSG00000128849'];assert len(target)==1
v={s:float(np.log2(float(target[s].iloc[0])/float(raw[s].sum())*1e6+.5)) for s in meta.index}
base=pd.read_csv(O.parent/'results/bulk_patient_deltas.tsv',sep='\t');expected=pd.read_csv(O/'results/GSE314812_GC05_CPM_sensitivity.tsv',sep='\t').set_index('variant');qa=[]
for variant in expected.index:
 yy=[];zz=[]
 for row in base.itertuples():
  if row.patient=='GC05' and variant=='exclude_GC05':continue
  pairs=[(row.PT_sample,row.PM_sample)]
  if row.patient=='GC05':
   if variant=='replace_with_b':pairs=[('GC05b_GCPT','GC05b_GCPM')]
   if variant=='average_a_b_within_patient':pairs=[('GC05a_GCPT','GC05a_GCPM'),('GC05b_GCPT','GC05b_GCPM')]
  yy.append(np.mean([v[b]-v[a] for a,b in pairs]));zz.append(np.mean([int(meta.loc[b,'preservation method']=='FFPE')-int(meta.loc[a,'preservation method']=='FFPE') for a,b in pairs]))
 y=np.array(yy);x=np.column_stack([np.ones(len(y)),zz]);inv=np.linalg.inv(x.T@x);beta=inv@x.T@y;res=y-x@beta;h=np.einsum('ij,jk,ik->i',x,inv,x);meat=x.T@np.diag((res/(1-h))**2)@x;cov=inv@meat@inv;p=2*t.sf(abs(beta[0]/np.sqrt(cov[0,0])),len(y)-np.linalg.matrix_rank(x))
 # SciPy exact distribution is valid here only after confirming no zero or tied absolute deltas.
 assert np.all(y!=0) and len(set(abs(y)))==len(y)
 vals={'n_patients_assumed':len(y),'mean_delta':y.mean(),'positive':int((y>0).sum()),'negative':int((y<0).sum()),'paired_rank_P':wilcoxon(y,method='exact').pvalue,'preservation_adjusted_effect':beta[0],'preservation_adjusted_P':p}
 errors={k:abs(float(val)-float(expected.loc[variant,k])) for k,val in vals.items()};assert max(errors.values())<1e-10;qa.append({'variant':variant,'max_error':max(errors.values())})
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'checks':qa,'independent_methods':'raw count recalculation; SciPy exact signed rank; direct HC3 matrix calculation; original primary TMM results retained'}
(O/'results/GC05_NUMERICAL_QA_110.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
for folder in [O/'raw/GSE308231_ranges',O/'raw/verified_ranges']:
 files=list(folder.rglob('*.bin'));partials=list(folder.rglob('*.partial'));print('TRANSFER',folder.name,'verified_chunks',len(files),'verified_bytes',sum(f.stat().st_size for f in files),'partial_bytes',sum(f.stat().st_size for f in partials))
print('INDEX', (O/'metadata/INDEX_COMPLETE.json').read_text())
print('LEDGER', (O/'MECHANISM_EVIDENCE_LEDGER_v6.json').read_text()[:18000])
