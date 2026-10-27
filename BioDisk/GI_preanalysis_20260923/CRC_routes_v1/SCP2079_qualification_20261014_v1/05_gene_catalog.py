from pathlib import Path
import json,requests
R=Path(__file__).resolve().parent
j=json.loads((R/'sources/explore.txt').read_text(encoding='utf-8'))
print({k:v for k,v in j.items() if k not in ['uniqueGenes','annotationList']})
print('uniqueSymbols',len(j['uniqueGenes']))
spec=json.loads(Path(r'C:\Users\31937\CGNL1_new_cohorts_20261013_v1\delivery\SPECIFICATION.json').read_text(encoding='utf-8'))
genes=sorted(set(['CGNL1']+spec['ordinary_markers']+sum(spec['panels'].values(),[])+sum(spec['outcomes'].values(),[])))
print('Missing fixed panel genes:',sorted(set(genes)-set(j['uniqueGenes'])))
(R/'GENE_AVAILABILITY.json').write_text(json.dumps({'source':'public visualization gene catalog; not a matrix integrity check','genes':genes,'missing':sorted(set(genes)-set(j['uniqueGenes']))},indent=2),encoding='utf-8')