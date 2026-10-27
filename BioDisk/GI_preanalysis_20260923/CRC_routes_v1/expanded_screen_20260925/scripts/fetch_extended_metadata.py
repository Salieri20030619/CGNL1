from pathlib import Path
import urllib.request,json,concurrent.futures
A=Path(__file__).resolve().parent.parent;M=A/'metadata';op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
jobs=[('zenodo15234826.json','https://zenodo.org/api/records/15234826'),('zenodo17638058.json','https://zenodo.org/api/records/17638058'),('zenodo10046840.json','https://zenodo.org/api/records/10046840'),('mendeley.html','https://data.mendeley.com/datasets/yd7bb3shn5/1'),('mendeley_dataset.json','https://data.mendeley.com/public-api/datasets/yd7bb3shn5'),('stanford_http.html','http://dna-discovery.stanford.edu/research/datasets/')]
for acc in ['13655','13651','13652']:
 jobs.append(('E-MTAB-'+acc+'.sdrf.txt',f'https://www.ebi.ac.uk/biostudies/files/E-MTAB-{acc}/E-MTAB-{acc}.sdrf.txt'))
def f(j):
 name,u=j
 try:
  with op.open(urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0'}),timeout=40) as r:b=r.read(10000000)
  (M/name).write_bytes(b);print(name,len(b),flush=True);return {'file':name,'url':u,'bytes':len(b),'status':'ok'}
 except Exception as e:print(name,str(e),flush=True);return {'file':name,'url':u,'status':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(f,jobs))
(A/'EXTENDED_METADATA_FETCH.json').write_text(json.dumps(r,indent=2))
