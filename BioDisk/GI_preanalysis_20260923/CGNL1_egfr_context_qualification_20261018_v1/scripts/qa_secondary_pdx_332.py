from pathlib import Path
import json,csv,gzip
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/secondary_pdx_329';link=json.loads((b/'metadata/DONOR_LINKAGE_324.json').read_text());old=json.loads((o/'ANALYSIS_COMPLETE_329.json').read_text())['contrasts'];checks=[];errors=[];total=0
for phase in ['PDX_72h','PDX_6w']:
 z=np.load(o/(phase+'_input.npz'));x=z['counts'];genes=list(z['genes']);rs=json.loads((o/(phase+'_samples.json')).read_text())
 with gzip.open(b/'raw/GSE204805_merged_hs_mm.tsv.gz','rt',encoding='utf-8') as f:
  rd=csv.reader(f,delimiter='\t');h=next(rd);ix=[h.index(r['sample_code']) for r in rs];j=0
  for row in rd:
   if row[0].startswith('H_'):
    assert row[0]==genes[j] and np.array_equal(np.array([int(row[i]) for i in ix]),x[j]);j+=1
 total+=x.size;ref=(x.min(axis=1)>0)&(z['genes']!='H_CGNL1');lr=np.log(x[ref]);ratio=np.sort(np.exp(lr-lr.mean(axis=1)[:,None]),axis=0);n=len(ratio);sf=(ratio[(n-1)//2]+ratio[n//2])/2;sf/=np.exp(np.log(sf).mean());assert np.max(abs(sf-z['sizefactors']))<1e-10
 for method in ['median_ratio','human_CPM']:
  v=np.log2(x[genes.index('H_CGNL1')]/sf+1) if method=='median_ratio' else np.log2(x[genes.index('H_CGNL1')]/x.sum(axis=0)*1e6+1);donors={}
  for case,donor in link[phase]['case_to_donor'].items():
   if donor:
    d=np.mean([v[i] for i,r in enumerate(rs) if r['patient_key']==case and r['treatment']=='CET'])-np.mean([v[i] for i,r in enumerate(rs) if r['patient_key']==case and r['treatment']=='NT']);donors.setdefault(donor,[]).append(d)
  d=np.array([np.mean(v) for v in donors.values()]);rank=rankdata(abs(d[d!=0]));obs=abs(sum(rank*np.sign(d[d!=0])));sums=np.array([0.])
  for r in rank:sums=np.concatenate([sums-r,sums+r])
  p=float(np.mean(abs(sums)>=obs-1e-12));r=next(r for r in old if r['phase']==phase and r['method']==method and r['gene']=='CGNL1');assert abs(p-r['P'])<1e-12;errors.append(abs(d.mean()-r['known_donors_mean']));checks.append({'phase':phase,'method':method,'n_donors':len(d),'P':p})
assert max(errors)<1e-10
res={'all_pass':True,'raw_count_values_checked':int(total),'exact_enumeration_verified':checks,'max_effect_error':max(errors),'scope':'Known-source-donor subsets only; unmapped cases descriptive'};(o/'INDEPENDENT_QA_332.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
