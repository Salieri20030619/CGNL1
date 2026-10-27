from pathlib import Path
import os,signal,json,datetime,requests,concurrent.futures,gzip,hashlib
O=Path.cwd()/'mechanism_investigation_v2';stopped=[]
for p in Path('/proc').iterdir():
 if p.name.isdigit():
  try:
   args=(p/'cmdline').read_bytes().split(b'\0')
   if any(a.endswith(b'/cgnl_mech44_acquire_route.py') for a in args):
    record={'pid':int(p.name),'starttime':(p/'stat').read_text().split()[21],'script':'cgnl_mech44_acquire_route.py'};os.kill(int(p.name),signal.SIGTERM);stopped.append(record)
  except (FileNotFoundError,PermissionError):pass
(O/'metadata/FTP_ROUTE_INTERVENTION.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'multiple gzip failures; repeated samefile FTP downloads diverge midstream; no defective input accepted','stopped_owned_processes':stopped,'preserved':'all intact andfailedfiles and receipts'},indent=2),encoding='utf-8')
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223165/suppl/GSE223165_fpkm_sample.txt.gz'
urls={'https_query_download.gz':base+'?download=1','http_sample.gz':base.replace('https:','http:')}
def get(kv):
 n,u=kv;r={'name':n,'url':u}
 try:
  a=requests.get(u,timeout=(15,60),headers={'Accept-Encoding':'identity'});p=O/'raw'/n;p.write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest());a.raise_for_status();gzip.decompress(a.content);r['gzip_crc']=True
 except Exception as e:r['error']=repr(e)
 print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(get,urls.items()))
(O/'metadata/HTTPS_FALLBACK_50.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
