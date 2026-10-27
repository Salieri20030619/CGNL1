from pathlib import Path
import csv,gzip,json,datetime,zipfile,xml.etree.ElementTree as E,re
import numpy as np
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024')
OUT=P/'Set03_mapping'
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Normal1-Normal9 mapping only, no disease comparison','features':'reuse11312nonCGNL1qualifiedfeatures; require positive normal export','anchors':{'processed1':10,'processed2':12,'processed3':13},'split':'reuse frozen hash halves','gate':'all three anchors and both halves agree for each sample; n>=500 rho>=0.98 margin>=0.02; accepted mappings unique; leftover channels unknown','biological_pairing_not_inferred':True}
f=OUT/'NORMAL_MAPPING_SPEC_BEFORE_RESULTS.json'
if not f.exists():f.write_text(json.dumps(spec,indent=2),encoding='utf-8')
with gzip.open(OUT/'CHANNEL_FEATURES.csv.gz','rt') as h:features=list(csv.DictReader(h))
lookup={x['feature']:i for i,x in enumerate(features)};normal=np.full((len(features),9),np.nan)
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(P/'hep-84-396-s006.xlsx') as z:
 ss=[''.join(x.itertext()) for x in E.fromstring(z.read('xl/sharedStrings.xml'))];rows=E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:sheetData/s:row',ns)
 def vals(row):
  d={}
  for c in row:
   v=c.find('s:v',ns)
   if v is not None:d[re.sub('[0-9]','',c.attrib['r'])]=ss[int(v.text)] if c.attrib.get('t')=='s' else v.text
  return d
 header={v:k for k,v in vals(rows[1]).items()}
 for row in rows[2:]:
  x=vals(row);m=re.search(r' ([STY])(\d+)$',x.get('A',''))
  if not m:continue
  key='|'.join([x.get('C',''),m[1],m[2]])
  if key not in lookup:continue
  for n in range(9):
   try:v=float(x.get(header[f'Normal{n+1}'],'nan'))
   except ValueError:continue
   if v>0:normal[lookup[key],n]=v
raw=np.log2(np.array([[float(x[f'raw{i}']) for i in range(1,17)] for x in features]));proc=np.log2(np.array([[float(x[f'processed{i}']) for i in range(1,4)] for x in features]));split=np.array([int(x['split']) for x in features]);results=[]
for n in range(9):
 for a,channel in enumerate([10,12,13]):
  for half in [0,1]:
   mask=(split==half)&np.isfinite(normal[:,n]);vv=np.log2(normal[mask,n])-proc[mask,a];cs=[]
   for c in range(16):
    if c+1==channel:continue
    rho=float(spearmanr(raw[mask,c]-raw[mask,channel-1],vv).statistic);cs.append((rho,c+1))
   cs.sort(reverse=True);results.append({'normal':f'Normal{n+1}','anchor_channel':channel,'split':half,'features':int(mask.sum()),'best_channel':cs[0][1],'rho':cs[0][0],'runnerup_rho':cs[1][0],'margin':cs[0][0]-cs[1][0]})
mapping={};accepted={}
for n in range(9):
 name=f'Normal{n+1}';rr=[x for x in results if x['normal']==name];channels={x['best_channel'] for x in rr};ok=len(channels)==1 and all(x['features']>=500 and x['rho']>=.98 and x['margin']>=.02 for x in rr);accepted[name]=ok
 if ok:mapping[name]=channels.pop()
unique=len(set(mapping.values()))==len(mapping) and not(set(mapping.values())&{10,12,13})
out={'results':results,'accepted_each':accepted,'mapping':mapping,'unique':unique,'all_pass':all(accepted.values()) and unique,'unassigned_channels':sorted(set(range(1,17))-{10,12,13}-set(mapping.values())),'normal_to_patient_pairing_unknown':True,'mapping_is_numerically_inferred':True}
(OUT/'NORMAL_MAPPING_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
with gzip.open(OUT/'NORMAL_MAPPING_FEATURES.csv.gz','wt',encoding='utf-8',newline='') as h:
 w=csv.writer(h);w.writerow(['feature','split']+[f'raw{i}' for i in range(1,17)]+[f'processed{i}' for i in range(1,4)]+[f'normal{i}' for i in range(1,10)])
 for k,x in enumerate(features):w.writerow([x['feature'],int(split[k])]+list(raw[k])+list(proc[k])+[np.log2(v) if np.isfinite(v) else '' for v in normal[k]])
print(json.dumps({k:v for k,v in out.items() if k!='results'}));print('rho_min',min(x['rho'] for x in results),'margin_min',min(x['margin'] for x in results))
