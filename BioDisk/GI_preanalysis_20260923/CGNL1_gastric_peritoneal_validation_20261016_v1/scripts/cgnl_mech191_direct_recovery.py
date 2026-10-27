from pathlib import Path
import os,json,datetime,hashlib,requests,concurrent.futures,threading,time,shutil
O=Path.cwd()/'mechanism_investigation_v2';probe=json.loads((O/'metadata/ENA_PROXY_PROBE_183.json').read_text(encoding='utf-8'));assert probe['success'] and probe['identical_to_verified_file_prefix'];plan=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text(encoding='utf-8'));lookup={t['name']:t for t in plan['tasks']};snapshot=json.loads((O/'metadata/LOCAL_FAILED_FILE_SNAPSHOT_186.json').read_text(encoding='utf-8'));tasks=[]
for r in snapshot['files']:
 assert not r['success'];t=lookup[r['name']];assert t['bytes']==r['bytes'] and t['md5']==r['md5'];tasks.append(t)
tasks=[t for t in tasks if t['name']=='DRR400687_1.fastq.gz'];assert len(tasks)==1
tasks.sort(key=lambda t:({'combined':0,'TGF':1,'IL1':2,'control':3}[t['condition']],{'GC41':0,'GC40':1,'GC42':2}[t['donor']] if t['condition']=='combined' else {'GC42':0,'GC41':1,'GC40':2}[t['donor']],t['name']))
R=O/'raw/recovery_direct_191';R.mkdir(exist_ok=True);assert shutil.disk_usage(R).free>sum(t['bytes'] for t in tasks)*2+5_000_000_000;tls=threading.local();lock=threading.Lock();receipts=[]
def save(p,x):
 q=p.with_name(p.name+'.writing');q.write_text(json.dumps(x,indent=2),encoding='utf-8');q.replace(p)
def digest(p,kind='md5'):
 h=hashlib.new(kind)
 with p.open('rb') as f:
  while b:=f.read(4194304):h.update(b)
 return h.hexdigest()
def session(route='proxy'):
 if not hasattr(tls,route):
  obj=requests.Session();obj.trust_env=False
  if route=='proxy':pass # verified direct ENA HTTPS; no proxy dependency
  setattr(tls,route,obj)
 return getattr(tls,route)
save(O/'metadata/FAILED_FILE_RECOVERY_PLAN_191.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'tasks':tasks,'workers':8,'source':'officialENAHTTPS; verifiedownedproxy primary, directHTTPS fallback on secondattempt','prior_live_writers_unchanged':True,'destination':'separate recoveryfolder; wholeMD5 thenatomicmasterhardlink','retry':'bounded3attempts with60/120secondbackoff andRetry-After respect; retainallpartials'})
for process in Path('/proc').iterdir():
 if process.name.isdigit():
  try:
   assert str(O.parent/'scripts/cgnl_mech177_server_failed_file_recovery.py').encode() not in (process/'cmdline').read_bytes().split(b'\0'),'Old177 writer still live'
  except (FileNotFoundError,PermissionError):pass
for t in tasks:
 rec=dict(t,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());master=O/'raw/DRA014736'/t['name'];folder=R/(t['name']+'.ranges');folder.mkdir(exist_ok=True);manifest={}
 try:
  if master.exists():
   assert master.stat().st_size==t['bytes'] and digest(master)==t['md5'];rec.update(success=True,reused_verified_master=True)
  else:
   # Test route before queueing the file, avoiding a burst of failing chunk requests.
   for attempt in range(3):
    try:
     with session('direct' if attempt==1 else 'proxy').get(t['url'],headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},stream=True,timeout=(20,45)) as a:assert a.status_code==206 and a.headers.get('Content-Range')==f'bytes 0-0/{t["bytes"]}'
     break
    except Exception:
     if attempt==2:raise
     time.sleep(60*(attempt+1))
   def piece(iv):
    start,end=iv;p=folder/f'{start}_{end}.bin';old=O/'raw/verified_ranges'/t['name']/p.name;oldman=old.parent/'range_manifest.json'
    prior=O/'raw/recovery_resumed_185'/(t['name']+'.ranges')/p.name;prior_manifest=prior.parent/'manifest.json'
    if prior.exists() and prior_manifest.exists():
     m=json.loads(prior_manifest.read_text(encoding='utf-8')).get(str(start),{})
     if m.get('success') and prior.stat().st_size==end-start and digest(prior,'sha256')==m.get('sha256'):
      shutil.copy2(prior,p)
      with lock:manifest[str(start)]={'bytes':end-start,'sha256':m['sha256'],'success':True,'reused_177_chunk':True};save(folder/'manifest.json',manifest)
      return p
    if old.exists() and oldman.exists():
     m=json.loads(oldman.read_text(encoding='utf-8')).get(str(start),{})
     if m.get('success') and old.stat().st_size==end-start and digest(old,'sha256')==m.get('sha256'):
      shutil.copy2(old,p)
      with lock:manifest[str(start)]={'bytes':end-start,'sha256':m['sha256'],'success':True,'reused_old_chunk':True};save(folder/'manifest.json',manifest)
      return p
    for attempt in range(1,4):
     tmp=p.with_name(p.name+f'.attempt{attempt}.partial');retry_after=0
     try:
      begin=time.monotonic();n=0;h=hashlib.sha256()
      with session('direct' if attempt==2 else 'proxy').get(t['url'],headers={'Range':f'bytes={start}-{end-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,120)) as a:
       if a.headers.get('Retry-After','').isdigit():retry_after=int(a.headers['Retry-After'])
       assert a.status_code==206 and a.headers.get('Content-Range')==f'bytes {start}-{end-1}/{t["bytes"]}'
       with tmp.open('wb') as f:
        for b in a.iter_content(65536):
         f.write(b);n+=len(b);h.update(b)
         if n>end-start or time.monotonic()-begin>3600:raise TimeoutError('onehour chunk bound')
      assert n==end-start;tmp.rename(p)
      with lock:manifest[str(start)]={'bytes':n,'sha256':h.hexdigest(),'success':True,'attempt':attempt};save(folder/'manifest.json',manifest)
      return p
     except Exception as e:
      err=repr(e)
      if attempt<3:time.sleep(max(60*attempt,retry_after))
    with lock:manifest[str(start)]={'success':False,'error':err};save(folder/'manifest.json',manifest)
    raise RuntimeError(err)
   ivs=[(s,min(s+8388608,t['bytes'])) for s in range(0,t['bytes'],8388608)]
   with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:parts=list(ex.map(piece,ivs))
   full=R/t['name'];tmp=full.with_name(full.name+'.assembled.partial');h=hashlib.md5();n=0
   with tmp.open('wb') as f:
    for p in parts:
     with p.open('rb') as src:
      while b:=src.read(4194304):f.write(b);n+=len(b);h.update(b)
   assert n==t['bytes'] and h.hexdigest()==t['md5'];tmp.rename(full)
   try:os.link(full,master);created=True
   except FileExistsError:assert master.stat().st_size==t['bytes'] and digest(master)==t['md5'];created=False
   rec.update(success=True,actual_md5=h.hexdigest(),master_hardlink_created=created)
 except Exception as e:rec.update(success=False,error=repr(e))
 rec['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipts.append(rec);save(O/'metadata/FAILED_FILE_RECOVERY_RECEIPTS_191.json',receipts);print(t['name'],rec.get('success'),rec.get('error','officialMD5verified'),flush=True)
save(O/'metadata/FAILED_FILE_RECOVERY_COMPLETE_191.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':len(receipts)==len(tasks) and all(r['success'] for r in receipts),'files':len(receipts),'scope':'only individually failed local files, notallstudyinputs'})
