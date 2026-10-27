from pathlib import Path
import json,csv,gzip,zipfile,xml.etree.ElementTree as E,re,hashlib,datetime,itertools
import numpy as np
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024');O=P/'Set03_protein_mapping';O.mkdir(exist_ok=True)
samples=['iCCA27','iCCA29','iCCA30','Normal5','Normal6']
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'protein channel identity and common-scaling provenance;noCGNL1effect','samples':samples,'features':'unique canonical protein accession;noCGNL1,contaminant,reverse,site-only;positive16raw and5export columns','split':'sha256accession parity','mapping_gate':'all10pairmaps agree across halves; n>=500 rho>=.98 margin>=.02','scaling':'record median log2(export ratio/raw ratio) and residual absolute95percentile on heldout;diagnostic only,not approved normalization','new_biological_tests':0}
f=O/'SPEC_BEFORE_RESULTS.json'
if not f.exists():f.write_text(json.dumps(spec,indent=2),encoding='utf-8')
raw={}
with (P/'Set03_GP/proteinGroups.txt').open() as h:
 for x in csv.DictReader(h,delimiter='\t'):
  if x['Gene names']=='CGNL1' or 'Q0VF96' in x['Protein IDs'] or any(x[k] for k in ['Reverse','Potential contaminant','Only identified by site']):continue
  key=x['Protein IDs']
  if ';' in key:continue
  v=[float(x[f'Reporter intensity corrected {i}'] or 0) for i in range(1,17)]
  if min(v)>0:raw.setdefault(key,[]).append(v)
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};joined=[]
with zipfile.ZipFile(P/'hep-84-396-s005.xlsx') as z:
 ss=[''.join(x.itertext()) for x in E.fromstring(z.read('xl/sharedStrings.xml'))];rows=E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:sheetData/s:row',ns)
 def vals(row):
  d={}
  for c in row:
   v=c.find('s:v',ns)
   if v is not None:d[re.sub('[0-9]','',c.attrib['r'])]=ss[int(v.text)] if c.attrib.get('t')=='s' else v.text
  return d
 header={v:k for k,v in vals(rows[1]).items()}
 for row in rows[2:]:
  x=vals(row);key=x.get('B')
  if key not in raw or len(raw[key])!=1:continue
  try:v=[float(x.get(header[s],'nan')) for s in samples]
  except ValueError:continue
  if all(np.isfinite(v)) and min(v)>0:joined.append((key,raw[key][0],v))
splits=np.array([int(hashlib.sha256(x[0].encode()).hexdigest(),16)%2 for x in joined]);a=np.log2(np.array([x[1] for x in joined]));b=np.log2(np.array([x[2] for x in joined]));pairs=list(itertools.permutations(range(16),2));results=[]
with gzip.open(O/'FEATURES.csv.gz','wt',encoding='utf-8',newline='') as h:
 w=csv.writer(h);w.writerow(['feature','split']+[f'raw{i}' for i in range(1,17)]+samples)
 for x,s in zip(joined,splits):w.writerow([x[0],int(s)]+x[1]+x[2])
for s in [0,1]:
 m=splits==s;assert sum(m)>=500
 rd=np.array([a[m,i]-a[m,j] for i,j in pairs]);rr=rankdata(rd,axis=1);rr-=rr.mean(axis=1,keepdims=True);rr/=np.linalg.norm(rr,axis=1,keepdims=True)
 for i,j in itertools.combinations(range(5),2):
  v=rankdata(b[m,i]-b[m,j]);v-=v.mean();v/=np.linalg.norm(v);cc=rr@v;ix=np.argsort(cc)[::-1];best,runner=ix[:2]
  results.append({'split':s,'sample_pair':[samples[i],samples[j]],'best_raw_channels':[v+1 for v in pairs[best]],'rho':float(cc[best]),'margin':float(cc[best]-cc[runner])})
mapping={};consistent=True
for x in results:
 for sample,ch in zip(x['sample_pair'],x['best_raw_channels']):
  if sample in mapping and mapping[sample]!=ch:consistent=False
  mapping[sample]=ch
passed=consistent and len(set(mapping.values()))==5 and all(x['rho']>=.98 and x['margin']>=.02 for x in results)
scaling=[]
if passed:
 for j in range(1,5):
  residual=(b[:,j]-b[:,0])-(a[:,mapping[samples[j]]-1]-a[:,mapping[samples[0]]-1]);train=float(np.median(residual[splits==0]));held=residual[splits==1]-train
  scaling.append({'sample':samples[j],'anchor':samples[0],'training_log2_scaling_offset':train,'heldout_median_residual':float(np.median(held)),'heldout_abs_residual_p95':float(np.quantile(np.abs(held),.95))})
out={'features':len(joined),'results':results,'mapping':mapping if consistent else None,'accepted':passed,'scaling_diagnostics_only':scaling,'no_CGNL1_effect_computed':True}
(O/'RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k!='results'}))
