from pathlib import Path
import requests,json,hashlib,datetime,re
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();D=O/'literature'
sources={
'macro_micro_2026':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13390483/fullTextXML',
'rim_2024':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11443209/fullTextXML',
'NM2_2023':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10202830/fullTextXML',
'endothelial_NM2_2025':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12238879/fullTextXML',
'HGP_2026':'https://www.nature.com/articles/s41416-026-03567-y',
'anergy_2026':'https://link.springer.com/article/10.1186/s12943-026-02732-4',
'GSE294385':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385&targ=self&form=text&view=full',
'GSE217414':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217414&targ=self&form=text&view=full',
'GSE335552':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335552&targ=self&form=text&view=full',
}
def fetch(item):
 name,url=item;rec={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  s=requests.Session();s.trust_env=False;r=s.get(url,timeout=35);b=r.content;rec.update(status=r.status_code,final_url=r.url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest());(D/(name+'.response')).write_bytes(b)
  if r.status_code==200:
   text=r.text
   from bs4 import BeautifulSoup
   soup=BeautifulSoup(text,'xml' if name in ['macro_micro_2026','rim_2024','NM2_2023','endothelial_NM2_2025'] else 'html');plain=soup.get_text(' ',strip=True);(D/(name+'.txt')).write_text(plain,encoding='utf-8')
   patterns=['Data availability','data availability','Zenodo','zenodo','GSE','CGNL1','histopathological','github']
   excerpts=[]
   for pat in patterns:
    for m in list(re.finditer(re.escape(pat),plain))[:8]:excerpts.append(plain[max(0,m.start()-80):m.start()+350])
   rec['inspection_excerpts']=excerpts
 except Exception as e:rec['error']=str(e)
 return rec
with ThreadPoolExecutor(max_workers=4) as pool:receipts=list(pool.map(fetch,sources.items()))
(D/'FETCH_RECEIPTS.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2),encoding='utf-8')
for r in receipts:print(json.dumps(r,ensure_ascii=False))
