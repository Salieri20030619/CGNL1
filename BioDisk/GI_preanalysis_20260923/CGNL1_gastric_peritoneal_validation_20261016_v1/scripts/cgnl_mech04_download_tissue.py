from pathlib import Path
import requests,json,hashlib,time,datetime,tarfile
O=Path.cwd()/'mechanism_investigation_v2';p=O/'raw/GSE308231_RAW.tar';tmp=p.with_suffix('.tar.partial')
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308231/suppl/GSE308231_RAW.tar'
r={'url':u,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
try:
 with requests.get(u,stream=True,timeout=(20,120)) as a:
  a.raise_for_status();r['header_size']=a.headers.get('Content-Length');h=hashlib.sha256();size=0;last=time.monotonic()
  with tmp.open('wb') as f:
   for chunk in a.iter_content(4*1024*1024):
    if chunk:f.write(chunk);h.update(chunk);size+=len(chunk)
    if time.monotonic()-last>60:print('downloaded',size,flush=True);last=time.monotonic()
  if r['header_size']:assert size==int(r['header_size'])
 with tarfile.open(tmp) as t:members=[{'name':m.name,'size':m.size} for m in t.getmembers()]
 assert len(members)>0
 tmp.rename(p);r.update(success=True,bytes=size,sha256=h.hexdigest(),members=members)
except Exception as e:r.update(success=False,error=str(e))
r['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'metadata/GSE308231_DOWNLOAD.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True)
if not r['success']:raise RuntimeError(r['error'])
