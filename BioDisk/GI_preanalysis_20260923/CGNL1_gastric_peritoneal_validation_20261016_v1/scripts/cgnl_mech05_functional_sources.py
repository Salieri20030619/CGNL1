from pathlib import Path
import requests,json,concurrent.futures,hashlib,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2'
ids=['PMC2555940','PMC10184829','PMC10202830','PMC13163822','PMC11513637','PMC4072573','PMC12159144','PMC13206751']
def get(i):
 u=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{i}/fullTextXML';r={'pmcid':i,'url':u}
 try:
  a=requests.get(u,timeout=(15,50));p=O/'literature'/(i+'.xml');p.write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest());a.raise_for_status();root=ET.fromstring(a.content)
  r['title']=' '.join(root.find('.//article-title').itertext());r['abstracts']=[' '.join(e.itertext()) for e in root.findall('.//abstract')]
  r['availability']=[' '.join(e.itertext()) for e in root.iter('sec') if e.find('title') is not None and 'availability' in ''.join(e.find('title').itertext()).lower()]
  print(i,r['title'],r['availability'],flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(get,ids))
(O/'metadata/FUNCTIONAL_SOURCE_CATALOG.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
