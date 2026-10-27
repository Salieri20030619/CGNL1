from pathlib import Path
import json,re
O=Path.cwd()/'mechanism_investigation_v2';s=(O/'metadata/GSE237876_samples.soft').read_bytes().decode('utf-8');out=[]
for b in s.split('^SAMPLE = ')[1:]:
 ls=b.splitlines();title=next(l.split(' = ',1)[1] for l in ls if l.startswith('!Sample_title ='));files=[l.split(' = ',1)[1] for l in ls if l.startswith('!Sample_supplementary_file')];patient=re.search(r'pat\d+',title);out.append({'accession':ls[0],'title':title,'patient':patient.group() if patient else None,'files':files})
(O/'metadata/ROUTE_INPUT_PLAN_43.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('SAMPLES',len(out),'PATIENTS',len(set(x['patient'] for x in out)));print(json.dumps(out[:3],indent=2));print('PERITONEAL_DONORS',sorted(set(x['patient'] for x in out if 'perit' in x['title'])))
