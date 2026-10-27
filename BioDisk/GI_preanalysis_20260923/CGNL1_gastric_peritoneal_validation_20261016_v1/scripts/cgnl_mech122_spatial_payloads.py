from pathlib import Path
import json,requests,xml.etree.ElementTree as ET,concurrent.futures,io,hashlib
from bs4 import BeautifulSoup
import openpyxl
O=Path.cwd()/'mechanism_investigation_v2';s=BeautifulSoup((O/'metadata/PMGC2026_publisher_121.raw').read_bytes(),'html.parser');text=s.get_text(' ',strip=True);print('ARTICLE_TEXT',text[:1800]);print('SECTIONS',[(x.get('id'),x.get_text(' ',strip=True)[:2500]) for x in s.find_all('section') if any(k in str(x.get('id')).lower() for k in ['data','sec','availability'])][-8:])
links=json.loads((O/'metadata/SPATIAL_ACCESS_121.json').read_text())[0]['data_links'];jobs=[x for x in links if x['href'].endswith('.xlsx')]
def get(x):
 a=requests.get(x['href'],timeout=(15,60));a.raise_for_status();p=O/'literature'/x['href'].rsplit('/',1)[1];p.write_bytes(a.content);wb=openpyxl.load_workbook(io.BytesIO(a.content),read_only=True,data_only=True);out={'url':x['href'],'bytes':len(a.content),'sha256':hashlib.sha256(a.content).hexdigest(),'sheets':[]}
 for ws in wb:
  rows=list(ws.values);hits=[{'row':i+1,'values':r} for i,r in enumerate(rows) if any('CGNL1' in str(v) for v in r)];out['sheets'].append({'sheet':ws.title,'rows':len(rows),'columns':ws.max_column,'first_rows':rows[:3],'target_matches':hits})
 wb.close();return out
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(get,jobs))
(O/'metadata/PMGC2026_SUPPLEMENT_122.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str),encoding='utf-8');print('SUPPLEMENT',json.dumps(out,ensure_ascii=False,default=str)[:15000])
for label in ['GC022','GC023']:
 root=ET.fromstring((O/f'metadata/PIPS_{label}_ENA_121.raw').read_bytes());e=root.find('.//EXPERIMENT_REF');assert e is not None;acc=e.attrib['accession'];a=requests.get('https://www.ebi.ac.uk/ena/browser/api/xml/'+acc,timeout=(15,45));a.raise_for_status();p=O/'metadata'/(acc+'_experiment.xml');p.write_bytes(a.content);print(label,acc,a.text[:6500])
