
from pathlib import Path
import requests,json,datetime,hashlib,gzip,concurrent.futures,threading,time,subprocess,os
p=Path.cwd();m=p/'mechanism_investigation_v2';dest=m/'raw/GSE286057_recovery_209';dest.mkdir(exist_ok=True)
tls=threading.local();recovered={}
def save(f,v):
 tmp=f.with_name(f.name+'.writing');tmp.write_text(json.dumps(v,indent=2),encoding='utf-8');tmp.replace(f)
def session():
 if not hasattr(tls,'s'):tls.s=requests.Session();tls.s.trust_env=False
 return tls.s
def digest(f):
 h=hashlib.sha256()
 with f.open('rb') as r:
  while b:=r.read(1048576):h.update(b)
 return h.hexdigest()
def recover(t):
 name=Path(t['file']).name;master=Path(t['file'])
 if name in recovered:return
 assert not master.exists(),'Unexpected destination exists; preserve and review'
 url=t['url']
 with session().get(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},stream=True,timeout=(20,60)) as r:
  assert r.status_code==206
  header=r.headers.get('Content-Range','');assert header.startswith('bytes 0-0/');size=int(header.split('/')[-1])
 folder=dest/(name+'.ranges');folder.mkdir(exist_ok=True)
 def piece(pair):
  a,b=pair;part=folder/(str(a)+'_'+str(b)+'.bin')
  for attempt in range(1,4):
   tmp=part.with_name(part.name+'.attempt'+str(attempt)+'.partial')
   try:
    with session().get(url,headers={'Range':'bytes='+str(a)+'-'+str(b-1),'Accept-Encoding':'identity'},stream=True,timeout=(20,90)) as r:
     assert r.status_code==206 and r.headers.get('Content-Range')=='bytes '+str(a)+'-'+str(b-1)+'/'+str(size)
     count=0
     with tmp.open('xb') as out:
      for block in r.iter_content(65536):out.write(block);count+=len(block)
    assert count==b-a,(name,a,b,count)
    tmp.rename(part);return part
   except Exception:
    if attempt==3:raise
    time.sleep(30*attempt)
 blocks=[(a,min(a+1048576,size)) for a in range(0,size,1048576)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:parts=list(ex.map(piece,blocks))
 full=dest/name
 with full.open('xb') as out:
  for part in parts:
   with part.open('rb') as src:
    while block:=src.read(1048576):out.write(block)
 assert full.stat().st_size==size
 with gzip.open(full,'rb') as inp:
  while inp.read(1048576):pass
 os.link(full,master)
 recovered[name]={'url':url,'expected_bytes':size,'bytes':size,'sha256':digest(full),'gzip_crc_verified':True,'success':True}
 save(m/'metadata/GSE286057_RECOVERY_RECEIPTS_209.json',recovered);print('RECOVERED',name,size,flush=True)
initial=json.loads((m/'metadata/GSE286057_DOWNLOAD_RECEIPTS_195.json').read_text(encoding='utf-8'))
for t in initial:
 if not t['success']:recover(t)
# Do not touch files still owned by the original downloader.
old=json.loads((m/'metadata/GSE286057_LAUNCH_195.json').read_text(encoding='utf-8'))
if not Path(old['receipt']).exists():subprocess.run(['/usr/bin/python3',str(p/'scripts/cgnl_wait_pidfd_system.py'),str(old['wrapper_pid'])],check=True)
assert Path(old['receipt']).exists(),'Original downloader exited without receipt'
final=json.loads((m/'metadata/GSE286057_DOWNLOAD_RECEIPTS_195.json').read_text(encoding='utf-8'));assert len(final)==108
for t in final:
 if not t['success']:recover(t)
valid=[]
for t in final:
 f=Path(t['file']);name=f.name;expected=recovered[name]['bytes'] if name in recovered else t['bytes']
 assert f.exists() and f.stat().st_size==expected
 with gzip.open(f,'rb') as inp:
  while inp.read(1048576):pass
 valid.append({'file':name,'bytes':expected,'sha256':digest(f),'gzip_crc_verified':True,'recovered':name in recovered})
save(m/'metadata/GSE286057_RECOVERED_INPUT_MANIFEST_209.json',valid)
save(m/'metadata/GSE286057_DOWNLOAD_RECOVERED_COMPLETE_209.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'files':108,'recovered_files':list(recovered),'prior_failures_and_partials_preserved':True})
print('ALL_108_INPUTS_VERIFIED',flush=True)
