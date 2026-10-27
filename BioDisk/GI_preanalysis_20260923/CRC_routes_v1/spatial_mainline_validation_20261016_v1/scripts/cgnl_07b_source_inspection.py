from pathlib import Path
import json,re,datetime
from bs4 import BeautifulSoup
O=Path.cwd();D=O/'literature';rows=[]
for name in ['macro_micro_PMC','rim_2024','HGP_2026','anergy_2026_retry','NM2_2023','endothelial_NM2_2025']:
 p=D/(name+'.txt')
 if not p.exists():continue
 text=p.read_text(encoding='utf-8');info={'source':name,'CGNL1_mentions':len(re.findall(r'\bCGNL1\b',text,re.I)),'accessions':sorted(set(re.findall(r'GSE\d+|E-MTAB-\d+|zenodo\.\d+',text))),'excerpts':[]}
 for pat in ['Data and code availability','Data availability','CGNL1','zenodo.org','17559845','19043057']:
  for m in list(re.finditer(re.escape(pat),text,re.I))[:2]:info['excerpts'].append(text[max(0,m.start()-40):m.start()+750])
 rows.append(info)
 print(name,'CGNL1 mentions',info['CGNL1_mentions'],'accessions',info['accessions'])
(D/'PRIMARY_SOURCE_INSPECTION.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
p=D/'HGP_project.response'
if p.exists():
 text=p.read_text(encoding='utf-8')
 try:
  j=json.loads(text);(O/'metadata/HGP_PROJECT_LINKS.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print('JSON PROJECT',list(j))
 except json.JSONDecodeError:
  soup=BeautifulSoup(text,'html.parser');links=[x.get('src') for x in soup.find_all('script')];fragments=[text[max(0,m.start()-100):m.start()+200] for m in re.finditer(r'iss_panel_2|\.h5ad|\.json|\.tmap',text)];record={'content_type':'HTML app response, not JSON project','scripts':links,'fragments':fragments};(O/'metadata/HGP_PROJECT_HTML_INSPECTION.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record)[:3500])
