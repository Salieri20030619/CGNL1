from pathlib import Path
import json
O=Path.cwd()/'mechanism_investigation_v2';rows=[]
for p in sorted((O/'metadata').glob('*_gastric_candidate.soft')):
 s=p.read_bytes().decode('utf-8');d={'id':p.name.split('_')[0]}
 for k in ['title','summary','overall_design','pubmed_id','supplementary_file']:d[k]=[l.split(' = ',1)[1] for l in s.splitlines() if l.startswith('!Series_'+k+' = ')]
 rows.append(d)
 print(d['id'],' | '.join(d['title']))
 if any(w in ' '.join(d['title']).lower() for w in ['gastric','peritoneal','gastritis']):print('DESIGN',' '.join(d['overall_design'])[:1000],'SUMMARY',' '.join(d['summary'])[:1100])
(O/'metadata/GASTRIC_CANDIDATE_DIGEST_42.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
