from pathlib import Path
import json,gzip,csv,itertools
import numpy as np
from scipy.stats import rankdata
import openpyxl
b=Path(__file__).resolve().parents[1];o=b/'results/independent_pdo_326';a=np.load(o/'input_and_normalization.npz');rs=json.loads((o/'samples.json').read_text());x=a['counts'];names=list(a['genes']);samples=list(a['samples']);errors=[]
with gzip.open(b/'raw/GSE204805_merged_hs_mm.tsv.gz','rt',encoding='utf-8') as f:
 r=csv.reader(f,delimiter='\t');h=next(r);ix=[h.index(z) for z in samples];j=0
 for row in r:
  if row[0].startswith('H_'):
   assert row[0]==names[j] and all(int(row[k])==int(x[j,i]) for i,k in enumerate(ix));j+=1
assert j==len(names)
refs=[j for j,g in enumerate(names) if g!='H_CGNL1' and min(x[j])>0]
logx=np.log(x[refs]);lograt=logx-logx.mean(axis=1)[:,None]
# median(exp(log-ratio)) agrees with original arithmetic median of ratios, including even-sized reference
ratios=np.sort(np.exp(lograt),axis=0);n=len(refs);sf=(ratios[(n-1)//2]+ratios[n//2])/2;sf=sf/np.exp(np.log(sf).mean());assert np.max(abs(sf-a['sizefactors']))<1e-10
keys=sorted(set(r['patient_key'] for r in rs));w=openpyxl.load_workbook(b/'raw/41467_2024_51909_MOESM4_ESM.xlsx',read_only=True,data_only=True)
mp={str(r[0]).strip():str(r[1]).strip() for r in w.worksheets[0].iter_rows(min_row=4,values_only=True) if r[0]};assert len({mp[k] for k in keys})==16
source=list(csv.DictReader((o/'patient_deltas.tsv').open(encoding='utf-8'),delimiter='\t'));checks=[]
for method in ['median_ratio','human_CPM']:
 v=np.log2(x[names.index('H_CGNL1')]/sf+1) if method=='median_ratio' else np.log2(x[names.index('H_CGNL1')]/x.sum(axis=0)*1e6+1)
 d=[]
 for k in keys:
  z=sum(v[i] for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='CET')/sum(r['patient_key']==k and r['treatment']=='CET' for r in rs)-sum(v[i] for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='NT')/sum(r['patient_key']==k and r['treatment']=='NT' for r in rs)
  old=float(next(r['delta'] for r in source if r['gene']=='CGNL1' and r['method']==method and r['case']==k));errors.append(abs(z-old));d.append(z)
 d=np.array(d);ranks=rankdata(abs(d[d!=0]));obs=abs(sum(ranks*np.sign(d[d!=0])));extreme=sum(abs(sum(ranks*np.array(signs)))>=obs-1e-12 for signs in itertools.product([-1,1],repeat=len(ranks)));p=extreme/2**len(ranks)
 old=next(r for r in json.loads((o/'ANALYSIS_COMPLETE_326.json').read_text())['contrasts'] if r['gene']=='CGNL1' and r['method']==method);assert abs(p-old['P'])<1e-12;checks.append({'method':method,'P':p})
assert max(errors)<1e-10
res={'all_pass':True,'raw_human_count_values_checked':int(x.size),'confirmed_distinct_donors':16,'normalization_reference_genes':len(refs),'max_delta_error':max(errors),'exact_test_independent_enumeration':checks}
(o/'INDEPENDENT_QA_327.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
