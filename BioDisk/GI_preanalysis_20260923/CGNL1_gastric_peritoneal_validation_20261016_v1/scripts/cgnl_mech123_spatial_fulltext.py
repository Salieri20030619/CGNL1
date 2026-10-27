from pathlib import Path
import json,requests,subprocess,hashlib,concurrent.futures,xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
from urllib.parse import urljoin
O=Path.cwd()/'mechanism_investigation_v2';s=BeautifulSoup((O/'metadata/PMGC2026_publisher_121.raw').read_bytes(),'html.parser');links=[urljoin('https://www.nature.com',x['href']) for x in s.find_all('a',href=True) if x.get_text(' ',strip=True)=='Download PDF'];assert links;u=links[0]
jobs=[('PMGC2026_accepted.pdf',u),('PIPS_GC022_sample.xml','https://www.ebi.ac.uk/ena/browser/api/xml/SRS22298254'),('PIPS_GC023_sample.xml','https://www.ebi.ac.uk/ena/browser/api/xml/SRS22298258')]
def get(job):
 name,u=job;a=requests.get(u,timeout=(15,90));a.raise_for_status();p=O/'literature'/name;p.write_bytes(a.content);r={'name':name,'url':u,'bytes':len(a.content),'sha256':hashlib.sha256(a.content).hexdigest()}
 if name.endswith('.pdf'):
  assert a.content.startswith(b'%PDF');t=p.with_suffix('.txt');subprocess.run(['pdftotext','-layout',str(p),str(t)],check=True);text=t.read_text(encoding='utf-8');r['hits']=[]
  for i,page in enumerate(text.split('\f')):
   if any(term in page.lower() for term in ['data availability','egad','processed data','figshare','zenodo','visium']):r['hits'].append({'page':i+1,'text':page})
 else:r['sample_metadata']=' '.join(ET.fromstring(a.content).itertext())
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(get,jobs))
(O/'metadata/SPATIAL_FULLTEXT_123.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
for r in out:
 if 'hits' in r:print('PMGC_DATA',json.dumps([x for x in r['hits'] if 'data availability' in x['text'].lower()],ensure_ascii=False)[:10000]);print('OTHER_HIT_PAGES',[x['page'] for x in r['hits']])
 else:print(json.dumps(r,ensure_ascii=False))
