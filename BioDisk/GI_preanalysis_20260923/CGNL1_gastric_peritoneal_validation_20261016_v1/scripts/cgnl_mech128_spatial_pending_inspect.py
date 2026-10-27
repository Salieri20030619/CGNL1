from pathlib import Path
import json,xml.etree.ElementTree as ET,requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
O=Path.cwd()/'mechanism_investigation_v2';rows=[]
for name in ['PIPS_GC022_sample.xml','PIPS_GC023_sample.xml']:
 p=O/'literature'/name
 if p.exists():
  root=ET.fromstring(p.read_bytes());r={'file':name,'title':root.findtext('.//TITLE'),'attributes':[{x.tag:x.text for x in a} for a in root.findall('.//SAMPLE_ATTRIBUTE')],'links':[{x.tag:x.text for x in a.iter() if x.text and x.text.strip()} for a in root.findall('.//SAMPLE_LINK')]};rows.append(r)
s=BeautifulSoup((O/'metadata/PMGC2026_publisher_121.raw').read_bytes(),'html.parser');links=[urljoin('https://www.nature.com',x['href']) for x in s.find_all('a',href=True) if x.get_text(' ',strip=True)=='Download PDF'];u=links[0]
try:
 a=requests.head(u,timeout=(15,30),allow_redirects=True);r={'url':u,'status':a.status_code,'headers':dict(a.headers)}
except Exception as e:r={'url':u,'error':repr(e)}
out={'sample_metadata':rows,'pdf_header_probe':r};(O/'metadata/SPATIAL_PENDING_128.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False)[:12000])
