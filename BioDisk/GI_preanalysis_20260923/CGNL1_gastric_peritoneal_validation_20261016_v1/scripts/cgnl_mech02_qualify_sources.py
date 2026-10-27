from pathlib import Path
import json,re
from bs4 import BeautifulSoup
O=Path.cwd()/'mechanism_investigation_v2';out={}
for name in ['ascites_europepmc','cgnl_ko_europepmc']:
 p=O/'literature'/(name+'.response');soup=BeautifulSoup(p.read_bytes(),'xml');text=soup.get_text(' ',strip=True)
 (O/'literature'/(name+'.txt')).write_text(text,encoding='utf-8')
 sections=[]
 for sec in soup.find_all('sec'):
  title=sec.find('title',recursive=False)
  if title and any(t in title.get_text().lower() for t in ['data availability','availability of data','supplementary','rna sequencing','rna-seq']):sections.append(sec.get_text(' ',strip=True))
 links=[dict(e.attrs) for e in soup.find_all(['ext-link','media','supplementary-material'])]
 ids=sorted(set(re.findall(r'(?:GSE\d+|PRJ\w+\d+|HRA\d+|CRA\d+|E-MTAB-\d+|OMIX\d+)',text)))
 out[name]={'sections':sections,'links':links,'ids':ids}
 print(name,'IDS',ids,'SECTIONS',sections)
 print('SUPPLEMENT_LINKS',[e for e in links if 'href' in str(e) and any(t in str(e) for t in ['xlsx','csv','zip','txt','GSE','CRA','OMIX','HRA'])])
for name in ['GSE308231','GSE308231_files','spatial_search']:
 t=(O/'literature'/(name+'.response')).read_text(encoding='utf-8');print(name,t);out[name]=t
(O/'metadata/SOURCE_QUALIFICATION_02.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
