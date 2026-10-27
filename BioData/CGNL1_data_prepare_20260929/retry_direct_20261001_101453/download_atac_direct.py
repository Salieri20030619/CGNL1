from pathlib import Path
import os,json,time,hashlib,urllib.request,urllib.error,concurrent.futures,threading,shutil,fcntl,traceback
R=Path(__file__).resolve().parent
lockfile=(R/'DOWNLOAD.lock').open('w');fcntl.flock(lockfile,fcntl.LOCK_EX|fcntl.LOCK_NB)
J=json.loads((R/'ATAC_RAW_MANIFEST_BLOCKED.json').read_text())['jobs'];L=threading.Lock();records={};progress={}
DIRECT=urllib.request.build_opener(urllib.request.ProxyHandler({}))
if (R/'MANIFEST.json').exists():
 prior=json.loads((R/'MANIFEST.json').read_text());planned={j['url'].rsplit('/',1)[-1]:j for j in J}
 for x in prior:
  name=x['url'].rsplit('/',1)[-1];p=Path(x.get('path',''));j=planned.get(name)
  if j and x.get('status') in ['downloaded_verified','verified_existing'] and x.get('sha256') and x.get('md5')==j['md5'] and p.is_file() and p.stat().st_size==j['bytes']:records[name]=x
def save(name,obj):
 p=R/name;t=R/(name+'.tmp');t.write_text(json.dumps(obj,indent=2));t.replace(p)
def hashes(p):
 md=hashlib.md5();sh=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):md.update(b);sh.update(b)
 return md.hexdigest(),sh.hexdigest()
def one(j):
 name=j['url'].rsplit('/',1)[-1];p=R/'ATAC_RAW'/j['run']/name;p.parent.mkdir(parents=True,exist_ok=True);part=Path(str(p)+'.part');last=None
 if p.exists():
  md,sha=hashes(p)
  if p.stat().st_size==j['bytes'] and md==j['md5']:return dict(j,path=str(p),status='verified_existing',sha256=sha)
  raise ValueError('existing final file differs from official checksum: '+str(p))
 for attempt in range(1,6):
  try:
   offset=part.stat().st_size if part.exists() else 0
   if offset>j['bytes']:raise ValueError('partial size exceeds expected')
   if offset<j['bytes']:
    if shutil.disk_usage(R).free-(j['bytes']-offset)<50*1024**3:raise RuntimeError('50GiB disk reserve reached')
    headers={'User-Agent':'CGNL1-data-preparation/1.1'}
    if offset:headers['Range']=f'bytes={offset}-'
    op=DIRECT
    with op.open(urllib.request.Request(j['url'],headers=headers),timeout=90) as resp:
     if offset and resp.status==206:
      if not resp.headers.get('Content-Range','').startswith('bytes '+str(offset)+'-'):raise ValueError('unexpected resume Content-Range')
     elif offset:offset=0
     with part.open('ab' if offset else 'wb') as f:
      n=offset;previous=0
      while True:
       b=resp.read(1024**2)
       if not b:break
       f.write(b);n+=len(b)
       if n>j['bytes']:raise ValueError('response exceeds official file size')
       if time.time()-previous>30:
        with L:
         progress[name]={'bytes':n,'expected':j['bytes'],'attempt':attempt,'updated':time.time()};save('ACTIVE_TRANSFERS.json',progress)
        previous=time.time()
   if part.stat().st_size!=j['bytes']:raise ValueError('download incomplete')
   md,sha=hashes(part)
   if md!=j['md5']:
    part.rename(Path(str(part)+'.badmd5.'+str(int(time.time()))));raise ValueError('official MD5 mismatch; retained failed payload')
   part.replace(p)
   return dict(j,path=str(p),status='downloaded_verified',sha256=sha,attempt=attempt,finished=time.time())
  except Exception as e:
   last=str(e)
   if 'disk reserve' in last:break
   if attempt<5:time.sleep(min(20*attempt,90))
 return dict(j,path=str(p),status='failed',error=last,finished=time.time())
def safe(j):
 try:return one(j)
 except Exception as e:return dict(j,status='failed',error=str(e))
save('PLAN.json',{'dataset':'E-MTAB-13652 / PRJEB71608','files':len(J),'bytes':sum(j['bytes'] for j in J),'workers':3,'reserve_GiB':50,'jobs':J})
save('STATE.json',{'status':'running','pid':os.getpid(),'started':time.time(),'files':len(J),'verified':len(records),'recorded':len(records),'network_mode':'direct_only'})
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 futures={ex.submit(safe,j):j for j in J if j['url'].rsplit('/',1)[-1] not in records}
 for future in concurrent.futures.as_completed(futures):
  x=future.result();name=x['url'].rsplit('/',1)[-1]
  with L:
   records[name]=x;progress.pop(name,None);save('MANIFEST.json',list(records.values()));save('ACTIVE_TRANSFERS.json',progress)
   save('STATE.json',{'status':'running','pid':os.getpid(),'recorded':len(records),'verified':sum(z['status'] in ['downloaded_verified','verified_existing'] for z in records.values()),'failed':sum(z['status']=='failed' for z in records.values()),'updated':time.time()})
  print(json.dumps(x),flush=True)
code=0 if all(x['status'] in ['downloaded_verified','verified_existing'] for x in records.values()) and len(records)==len(J) else 1
save('COMPLETE.json',{'exit_code':code,'planned':len(J),'verified':sum(x['status'] in ['downloaded_verified','verified_existing'] for x in records.values()),'failed':sum(x['status']=='failed' for x in records.values()),'finished':time.time(),'official_md5_required':True})
save('STATE.json',{'status':'complete' if code==0 else 'needs_attention','exit_code':code,'updated':time.time()})
raise SystemExit(code)
