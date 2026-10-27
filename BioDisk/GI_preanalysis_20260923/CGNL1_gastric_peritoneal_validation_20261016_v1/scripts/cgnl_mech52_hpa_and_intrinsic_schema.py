from pathlib import Path
from bs4 import BeautifulSoup
import re,json
O=Path.cwd()/'mechanism_investigation_v2';s=(O/'metadata/HPA_CGNL1_gastric.html').read_bytes().decode('utf-8');soup=BeautifulSoup(s,'html.parser');t=soup.get_text(' ',strip=True);out={}
out['antibody_ids']=sorted(set(re.findall(r'HPA\d{6}|CAB\d{6}',s)));out['sections']=[]
for pat in ['Antibody staining','Staining','Patients','Reliability','High','Not detected']:
 for m in list(re.finditer(re.escape(pat),t))[:3]:out['sections'].append({'term':pat,'text':t[max(0,m.start()-120):m.start()+900]})
out['image_links']=[i.get('src') for i in soup.find_all('img') if any(z in str(i.get('src')) for z in ['staining','images.proteinatlas','dictionary'])][:20]
(O/'metadata/HPA_GASTRIC_ANNOTATION_52.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False)[:7500])
for n in ['GSE53276.soft','GSE53276_samples.soft']:
 p=O/'metadata'/n
 if p.exists():
  t=p.read_bytes().decode('utf-8');ls=[l for l in t.splitlines() if l.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Series_supplementary_file','!Sample_title','!Sample_characteristics','!Sample_description','!Sample_platform_id'))];print(n,'\n'.join(ls)[:4500])
