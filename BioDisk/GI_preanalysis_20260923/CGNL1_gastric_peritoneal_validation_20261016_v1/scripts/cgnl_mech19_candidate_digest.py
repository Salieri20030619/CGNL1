from pathlib import Path
import json
O=Path.cwd()/'mechanism_investigation_v2';rows=[]
for p in sorted((O/'metadata').glob('*_candidate.soft')):
 text=p.read_bytes().decode('utf-8');row={'file':p.name}
 for key in ['title','overall_design','supplementary_file']:
  row[key]=[l.split(' = ',1)[1] for l in text.splitlines() if l.startswith('!Series_'+key+' = ')]
 rows.append(row)
 print(p.name,' | '.join(row['title']), 'DESIGN', ' '.join(row['overall_design'])[:420])
(O/'metadata/CANDIDATE_DIGEST_partial.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
