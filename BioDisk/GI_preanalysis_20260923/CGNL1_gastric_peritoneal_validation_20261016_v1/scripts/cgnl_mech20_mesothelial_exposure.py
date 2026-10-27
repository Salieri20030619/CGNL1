from pathlib import Path
import requests,urllib.request,json,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';ids=['GSE223165','GSE120303','GSE121372'];tasks=[]
for i in ids:
 p=O/'metadata'/(i+'_candidate.soft');s=p.read_bytes().decode('utf-8');print(i,'\n'.join(l for l in s.splitlines() if l.startswith(('!Series_title','!Series_pubmed_id','!Series_supplementary_file'))),flush=True)
 tasks.append((i+'_samples.soft','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+i+'&targ=gsm&form=text&view=full'))
 for l in s.splitlines():
  if l.startswith('!Series_supplementary_file = '):
   u=l.split(' = ',1)[1]
   if not u.endswith('_RAW.tar'):tasks.append((u.rsplit('/',1)[1],u))
def fetch(item):
 n,u=item;r={'file':n,'url':u}
 try:
  if u.startswith('ftp:'):
   with urllib.request.urlopen(u,timeout=60) as a:b=a.read(60*1024*1024+1)
   assert len(b)<=60*1024*1024,'defer larger supplement to dedicated downloader'
  else:
   a=requests.get(u,timeout=(15,60));a.raise_for_status();b=a.content
  p=O/('metadata' if n.endswith('.soft') else 'raw')/n;p.write_bytes(b);r.update(success=True,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),magic_hex=b[:8].hex());print(n,len(b),flush=True)
 except Exception as e:r.update(success=False,error=str(e));print(n,str(e),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(fetch,tasks))
(O/'metadata/ACQUISITION_20.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
