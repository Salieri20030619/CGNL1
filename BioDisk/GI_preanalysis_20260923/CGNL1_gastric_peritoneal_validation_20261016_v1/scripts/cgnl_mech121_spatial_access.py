from pathlib import Path
import requests,json,datetime,concurrent.futures,xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
O=Path.cwd()/'mechanism_investigation_v2'
jobs=[('PMGC2026_publisher','https://www.nature.com/articles/s42003-026-10571-8'),('DSP2024_publisher','https://els-jbs-prod-cdn.jbs.elsevierhealth.com/article/S0016-5085%2824%2905348-4/fulltext'),('PIPS_GC022_ENA','https://www.ebi.ac.uk/ena/browser/api/xml/SRR30190686'),('PIPS_GC023_ENA','https://www.ebi.ac.uk/ena/browser/api/xml/SRR30190685')]
def get(pair):
 name,u=pair;r={'name':name,'url':u}
 try:
  a=requests.get(u,timeout=(15,45));r.update(status=a.status_code,final_url=a.url);a.raise_for_status();p=O/'metadata'/(name+'_121.raw');p.write_bytes(a.content);s=BeautifulSoup(a.content,'html.parser');text=s.get_text(' ',strip=True);r['bytes']=len(a.content);r['hits']=[]
  for term in ['Data availability','Data Availability','data availability','EGAD','GSE','figshare','Visium','visium','CytAssist']:
   pos=text.find(term)
   if pos>=0:r['hits'].append({'term':term,'text':text[max(0,pos-200):pos+2200]})
  r['data_links']=[{'text':x.get_text(' ',strip=True),'href':x.get('href')} for x in s.find_all('a') if any(k in (str(x.get('href'))+' '+x.get_text()).lower() for k in ['supp','source data','figshare','zenodo','github','ega-archive','moesm'])]
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(get,jobs))
(O/'metadata/SPATIAL_ACCESS_121.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
for r in results:print(json.dumps(r,ensure_ascii=False)[:13000])
q=json.loads((O/'metadata/2025_PIPS_SPATIAL_QUALIFICATION_120.json').read_text(encoding='utf-8'));print('PIPS_DESIGN',*[s for s in q['paragraphs'] if s.startswith('Spatial transcriptomics were')],sep='\n')
