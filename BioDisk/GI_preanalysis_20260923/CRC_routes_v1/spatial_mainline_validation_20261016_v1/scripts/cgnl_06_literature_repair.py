from pathlib import Path
import requests,json,re,hashlib,datetime,xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();D=O/'literature';receipts=[]
for name in ['rim_2024','NM2_2023','endothelial_NM2_2025']:
 p=D/(name+'.response');root=ET.fromstring(p.read_bytes());plain=' '.join(root.itertext());(D/(name+'.txt')).write_text(plain,encoding='utf-8');receipts.append({'name':name,'repair':'standard-library XML parser; original response preserved','chars':len(plain)})
sources={'macro_micro_PMC':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/','zenodo_HD_current':'https://zenodo.org/api/records/19043057/versions/latest','zenodo_HD_initial':'https://zenodo.org/api/records/17559845','HGP_project':'https://crclm2.serve.scilifelab.se/iss_panel_2.tmap','HGP_code_tree':'https://api.github.com/repos/axanderssonumu/iss_crclm/git/trees/main?recursive=1','anergy_2026_retry':'https://link.springer.com/article/10.1186/s12943-026-02732-4'}
def fetch(item):
 name,url=item;r={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  s=requests.Session();s.trust_env=False;q=s.get(url,timeout=35);b=q.content;(D/(name+'.response')).write_bytes(b);r.update(status=q.status_code,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
  if q.status_code==200:
   if 'zenodo' in name:
    j=q.json();r.update(access=j.get('access',j.get('metadata',{}).get('access_right')),files=[{'key':a.get('key'),'size':a.get('size')} for a in j.get('files',[])],id=j.get('id'))
   elif name=='HGP_project':r['keys']=list(q.json());r['entries']=str(q.json())[:2500]
   else:
    text=BeautifulSoup(q.text,'html.parser').get_text(' ',strip=True);(D/(name+'.txt')).write_text(text,encoding='utf-8');r['CGNL1_mentions']=len(re.findall('CGNL1',text,re.I));r['accession_ids']=sorted(set(re.findall(r'GSE\d+|E-MTAB-\d+|zenodo\.\d+',text)))
    r['links']=[a.get('href') for a in BeautifulSoup(q.text,'html.parser').find_all('a') if any(k in str(a.get('href')) for k in ['zenodo','github','MOESM','mmc','geo/query'])]
 except Exception as e:r['error']=str(e)
 return r
with ThreadPoolExecutor(max_workers=4) as pool:receipts.extend(pool.map(fetch,sources.items()))
(D/'FETCH_REPAIR_RECEIPTS.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2),encoding='utf-8')
for r in receipts:print(json.dumps(r,ensure_ascii=False))
