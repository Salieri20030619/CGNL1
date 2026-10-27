from pathlib import Path
import subprocess,hashlib,json,concurrent.futures,threading,datetime,re,shlex,time
P=Path(r'C:\Users\31937\cgnl_transfer\restart_20261024');D=P/'C04_full';D.mkdir(exist_ok=True)
spec=json.loads((P/'C04_full_ACQUISITION_SPEC.json').read_text());files=[x for x in spec['files'] if x['name']!='phospho_fraction_1.raw']
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1';root=remote+'/continuous_exploration_20261023_v1/C04_full_fraction_qualification';script=remote+'/scripts/cgnl_c04_accept_convert_file.py'
lock=threading.Lock();state={'started':datetime.datetime.now(datetime.timezone.utc).isoformat(),'download_workers':1,'conversion_workers':2,'files':{},'reused':['phospho_fraction_1.raw'],'complete':False}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as stream:
  for b in iter(lambda:stream.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def update(name,**kw):
 with lock:
  state['files'].setdefault(name,{}).update(kw);tmp=D/'STATUS.tmp';tmp.write_text(json.dumps(state,indent=2),encoding='utf-8');tmp.replace(D/'STATUS.json')
 print(json.dumps(dict(file=name,**kw)),flush=True)
def run(cmd):return subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
def download(item):
 name=item['name'];assert re.fullmatch(r'phospho_(fraction_\d+|pY)\.raw',name);f=D/name;part=D/(name+'.partial');receipt=D/(name+'.receipt.json')
 if f.exists():
  assert receipt.exists();r=json.loads(receipt.read_text());assert f.stat().st_size==item['bytes'] and sha(f)==r['sha256'];return item,r
 attempts=[];update(name,stage='downloading')
 for attempt in range(2):
  cmd=['curl.exe','-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time','3600','--speed-limit','1024','--speed-time','180','--silent','--show-error']
  if part.exists() and part.stat().st_size:cmd+=['--continue-at','-']
  cmd+=['--output',str(part),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name];r=run(cmd);attempts.append(dict(exit=r.returncode,bytes=part.stat().st_size if part.exists() else 0,stderr=r.stderr))
  if r.returncode==0 and part.stat().st_size==item['bytes']:break
 else:
  update(name,stage='download_failed',attempts=attempts);return None
 with part.open('rb') as h:assert h.read(2)==b'\x01\xa1'
 digest=sha(part);part.rename(f);r=dict(name=name,bytes=f.stat().st_size,expected=item['bytes'],sha256=digest,route='Windowsnativecurl FTPS certificatevalidation',attempts=attempts);receipt.write_text(json.dumps(r,indent=2),encoding='utf-8');update(name,stage='download_verified',bytes=r['bytes'],sha256=digest);return item,r
def transfer_convert(pair):
 if pair is None:return False
 item,r=pair;name=item['name'];f=D/name;update(name,stage='transferring')
 cmd=['scp','-q','-o','BatchMode=yes','-o','ConnectTimeout=20',str(f),'sgs10:'+root+'/raw/'+name+'.uploading'];t=run(cmd)
 if t.returncode: update(name,stage='transfer_failed',exit=t.returncode,stderr=t.stderr);return False
 update(name,stage='converting')
 command='python3 '+shlex.quote(script)+' '+shlex.quote(name)+' '+shlex.quote(r['sha256']);t=run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10',command]);(D/(name+'.remote.log')).write_text(t.stdout+'\n'+t.stderr,encoding='utf-8')
 if t.returncode: update(name,stage='remote_conversion_or_QA_failed',exit=t.returncode,stderr=t.stderr[-800:]);return False
 z=json.loads(t.stdout.strip().splitlines()[-1]);assert z['raw_sha256']==r['sha256'];(D/(name+'.remote_acceptance.json')).write_text(json.dumps(z,indent=2),encoding='utf-8');update(name,stage='accepted' if z['acquisition_pass'] else 'acquisition_failed',spectra=z['spectrum_counts'],warnings=len(z['conversion_warnings']));return z['acquisition_pass']
pending=[];errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as downloads,concurrent.futures.ThreadPoolExecutor(max_workers=2) as converters:
 tasks={downloads.submit(download,f):f['name'] for f in files}
 for task in concurrent.futures.as_completed(tasks):
  try:
   pair=task.result()
   if pair is not None:pending.append(converters.submit(transfer_convert,pair))
   else:errors.append(tasks[task])
  except Exception as e:errors.append(tasks[task]);update(tasks[task],stage='exception',error=repr(e))
 for task in concurrent.futures.as_completed(pending):
  try:
   if not task.result():errors.append('transfer_or_conversion')
  except Exception as e:errors.append(repr(e));print(repr(e),flush=True)
state['complete']=not errors and len(state['files'])==12 and all(x['stage']=='accepted' for x in state['files'].values());state['finished']=datetime.datetime.now(datetime.timezone.utc).isoformat();state['errors']=errors;(D/'STATUS.json').write_text(json.dumps(state,indent=2),encoding='utf-8');print(json.dumps({'complete':state['complete'],'errors':errors}),flush=True)
