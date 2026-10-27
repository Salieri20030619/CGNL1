from pathlib import Path
import json,csv,gzip,zipfile,xml.etree.ElementTree as ET,re,hashlib,datetime,itertools
import numpy as np
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024')
OUT=P/'Set03_mapping';OUT.mkdir(exist_ok=True)
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'numerical provenance mapping, not biological effect testing','set':3,'exclude_gene':'CGNL1','features':'exact unique accession+residue+position; raw localization>=0.75; mono-phospho only; no contaminant/decoy; positive values in all16rawchannels and3knownset3iCCAtumors','split':'sha256(featurekey) parity, deterministic discovery and heldout halves','acceptance':'each of3oriented tumor-pair mappings consistent, each half n>=500, Spearman>=0.98 and margin>=0.02 over next orientedchannelpair in both halves','mapping_is_inferred_not_author_documented':True,'new_biological_tests':0}
f=OUT/'CHANNEL_AUDIT_SPEC_BEFORE_RESULTS.json'
if not f.exists():f.write_text(json.dumps(spec,indent=2),encoding='utf-8')
clinical=json.loads((P/'Mayo_CLINICAL_ROWS.json').read_text());samples=[x['B'] for x in clinical if x.get('G')=='3' and x.get('D')=='ICCA'];assert len(samples)==3
raw={}
with (P/'Set03_Phospho/Phospho (STY)Sites.txt').open() as h:
 for x in csv.DictReader(h,delimiter='\t'):
  if 'CGNL1' in x['Gene names'].split(';') or x['Reverse'] or x['Potential contaminant']:continue
  if x['Number of Phospho (STY)']!='1' or float(x['Localization prob'])<.75:continue
  if ';' in x['Proteins'] or not x['Positions within proteins'].isdigit():continue
  key=(x['Proteins'],x['Amino acid'],x['Positions within proteins'])
  vals=[float(x[f'Reporter intensity corrected {i}___1'] or 0) for i in range(1,17)]
  if min(vals)>0:raw.setdefault(key,[]).append(vals)
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
book=P/'hep-84-396-s006.xlsx';joined=[]
with zipfile.ZipFile(book) as z:
 strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
 root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'));rows=root.findall('.//s:sheetData/s:row',ns)
 def values(row):
  d={}
  for c in row:
   v=c.find('s:v',ns)
   if v is not None:d[re.sub('[0-9]','',c.attrib['r'])]=strings[int(v.text)] if c.attrib.get('t')=='s' else v.text
  return d
 header=values(rows[1]);columns={v:k for k,v in header.items()}
 for row in rows[2:]:
  x=values(row);match=re.search(r' ([STY])(\d+)$',x.get('A',''))
  if not match or x.get('B')=='CGNL1':continue
  key=(x.get('C'),match[1],match[2])
  if key not in raw or len(raw[key])!=1:continue
  try:v=[float(x.get(columns[s],'nan')) for s in samples]
  except ValueError:continue
  if all(np.isfinite(v)) and min(v)>0:joined.append((key,raw[key][0],v))
assert len(joined)>=1000,len(joined)
splits=np.array([int(hashlib.sha256('|'.join(x[0]).encode()).hexdigest(),16)%2 for x in joined])
with gzip.open(OUT/'CHANNEL_FEATURES.csv.gz','wt',encoding='utf-8',newline='') as h:
 w=csv.writer(h);w.writerow(['feature','split']+[f'raw{i}' for i in range(1,17)]+[f'processed{i}' for i in range(1,len(samples)+1)])
 for x,split in zip(joined,splits):w.writerow(['|'.join(x[0]),int(split)]+x[1]+x[2])
a=np.log2(np.array([x[1] for x in joined]));b=np.log2(np.array([x[2] for x in joined]));pairs=list(itertools.permutations(range(16),2));results=[]
for split in [0,1]:
 m=splits==split;assert sum(m)>=500
 rawdiff=np.array([a[m,i]-a[m,j] for i,j in pairs]);rr=rankdata(rawdiff,axis=1);rr-=rr.mean(axis=1,keepdims=True);rr/=np.linalg.norm(rr,axis=1,keepdims=True)
 for i,j in itertools.combinations(range(len(samples)),2):
  v=rankdata(b[m,i]-b[m,j]);v-=v.mean();v/=np.linalg.norm(v);cors=rr@v;order=np.argsort(cors)[::-1];best,second=order[:2];c1,c2=pairs[best]
  results.append({'split':split,'features':int(sum(m)),'sample_pair':[samples[i],samples[j]],'best_raw_channels':[c1+1,c2+1],'rho':float(cors[best]),'runnerup_rho':float(cors[second]),'margin':float(cors[best]-cors[second])})
mapping={};consistent=True
for x in results:
 for sample,channel in zip(x['sample_pair'],x['best_raw_channels']):
  if sample in mapping and mapping[sample]!=channel:consistent=False
  mapping[sample]=channel
passed=consistent and len(set(mapping.values()))==len(samples) and all(x['rho']>=.98 and x['margin']>=.02 for x in results)
out={'features_total':len(joined),'samples':samples,'results':results,'consistent':consistent,'accepted':passed,'provisional_mapping':mapping if consistent else None,'mapping_is_numerically_inferred':True,'no_biological_tests':True}
(OUT/'CHANNEL_AUDIT_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(out),flush=True)
