from pathlib import Path
import json,collections,gzip
O=Path.cwd()/'mechanism_investigation_v2';out={}
for i in ['GSE223165','GSE120303','GSE121372']:
 p=O/'metadata'/(i+'_samples.soft')
 if not p.exists():continue
 s=p.read_bytes().decode('utf-8');samples=[]
 for block in s.split('^SAMPLE = ')[1:]:
  ls=block.splitlines();d={'accession':ls[0]}
  for k in ['title','source_name_ch1','characteristics_ch1','treatment_protocol_ch1','growth_protocol_ch1','data_processing','platform_id']:
   d[k]=[l.split(' = ',1)[1] for l in ls if l.startswith('!Sample_'+k+' = ')]
  samples.append(d)
 out[i]=samples;print(i,'SAMPLES',len(samples))
 for d in samples:print(d['accession'],d['title'],d['characteristics_ch1'])
 if samples:print('PROTOCOL',samples[0])
for p in sorted((O/'raw').glob('GSE223165*.gz')):
 try:
  with gzip.open(p,'rt',encoding='utf-8') as f:head=f.readline();first=f.readline()
  print(p.name,'HEADER',head[:1800],'FIRST',first[:300])
 except Exception as e:print(p.name,str(e))
(O/'metadata/EXPOSURE_DESIGN_22.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
