from pathlib import Path
import gzip,json
O=Path.cwd()/'mechanism_investigation_v2';files=sorted((O/'raw/GSE237876').glob('*.txt.gz'));out=[]
for p in files[:3]:
 with gzip.open(p,'rt',encoding='utf-8') as f:
  header=f.readline();first=[f.readline().rstrip() for _ in range(3)]
 out.append({'file':p.name,'header':header,'first_rows':first})
s=(O/'metadata/GSE237876_samples.soft').read_bytes().decode('utf-8');methods=sorted(set(l for l in s.splitlines() if l.startswith('!Sample_data_processing')))
r={'complete_files_current':len(files),'schemas':out,'author_methods':methods};(O/'metadata/ROUTE_SCHEMA_69.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,ensure_ascii=False))
