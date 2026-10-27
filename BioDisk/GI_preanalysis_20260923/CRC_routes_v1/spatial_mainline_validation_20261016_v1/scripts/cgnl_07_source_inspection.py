from pathlib import Path
import json,re,datetime
O=Path.cwd();D=O/'literature';rows=[]
for name in ['macro_micro_PMC','rim_2024','HGP_2026','anergy_2026_retry']:
 p=D/(name+'.txt')
 if not p.exists():continue
 text=p.read_text(encoding='utf-8');info={'source':name,'CGNL1_mentions':len(re.findall(r'\bCGNL1\b',text,re.I)),'accessions':sorted(set(re.findall(r'GSE\d+|E-MTAB-\d+|zenodo\.\d+',text))),'excerpts':[]}
 for pat in ['Data and code availability','Data availability','CGNL1','zenodo.org','github.com','17559845','19043057']:
  for m in list(re.finditer(re.escape(pat),text,re.I))[:2]:info['excerpts'].append(text[max(0,m.start()-40):m.start()+750])
 rows.append(info);print(json.dumps(info,ensure_ascii=False))
p=D/'HGP_project.response'
if p.exists():
 j=json.loads(p.read_text(encoding='utf-8'));print('HGP_PROJECT_KEYS',list(j));print('HGP_PROJECT_FILES',json.dumps({k:v for k,v in j.items() if k in ['layers','markerFiles','regionFiles','collectionMode','filename']},ensure_ascii=False)[:6000]);(O/'metadata/HGP_PROJECT_LINKS.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'PRIMARY_SOURCE_INSPECTION.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
