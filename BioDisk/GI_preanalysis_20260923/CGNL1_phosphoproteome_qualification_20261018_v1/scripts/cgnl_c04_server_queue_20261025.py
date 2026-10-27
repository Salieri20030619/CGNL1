from pathlib import Path
import json,os,hashlib,subprocess,datetime,math,shutil
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
F=P/'continuous_exploration_20261023_v1/C04_full_fraction_qualification';W=F/'server_queue_20261025_v1'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,obj):
 t=path.with_name(path.name+'.tmp');t.write_text(json.dumps(obj,indent=2),encoding='utf-8');t.replace(path)
def sha(path,limit=None):
 h=hashlib.sha256();remaining=limit
 with path.open('rb') as f:
  while remaining is None or remaining>0:
   b=f.read(8*1024*1024 if remaining is None else min(remaining,8*1024*1024))
   if not b:break
   h.update(b)
   if remaining is not None:remaining-=len(b)
 assert remaining in (None,0);return h.hexdigest()
prefixes={'phospho_fraction_12.raw':(78675968,'b7bd1e35a7ff21a71794226bd9805edd50044ebbaa37e64cfbe6a5a5efa0c67b'),'phospho_fraction_2.raw':(6029312,'0b14d658438fc7849d2fc900cb5025dbb90e3d540ca80197eb0ca3a735226fe2')}
spec=json.loads((F/'ACQUISITION_SPEC_BEFORE_RESULTS.json').read_text());assert len(spec['files'])==13
W.mkdir(exist_ok=True);lock=W/'ACTIVE.json'
with lock.open('x',encoding='utf-8') as h:json.dump({'pid':os.getpid(),'utc':now(),'script_sha256':sha(Path(__file__))},h,indent=2)
state={'started':now(),'pid':os.getpid(),'files':{},'download_workers':1,'processing_workers':1,'whole_project_complete':False,'complete':False}
def update(name,**kw):
 state['files'].setdefault(name,{}).update(kw);save(W/'STATUS.json',state);print(json.dumps({'file':name,**kw}),flush=True)
try:
 assert shutil.disk_usage(F).free>15*1024**3,'Insufficient free disk for frozen scope'
 for name,(size,digest) in prefixes.items():
  f=F/'raw'/(name+'.server.partial');assert f.stat().st_size==size and sha(f)==digest
 save(W/'PREFIX_ACCEPTANCE.json',{n:{'bytes':s,'sha256':h,'matches_windows':True} for n,(s,h) in prefixes.items()})
 for item in spec['files']:
  name=item['name'];stem=Path(name).stem
  if name=='phospho_fraction_1.raw':
   assert (F/'REUSED_FIRST_FRACTION_SEARCH.json').exists();update(name,stage='reused_pilot');continue
  complete=F/'searches'/stem/'FILE_SEARCH_COMPLETE.json'
  if complete.exists():
   ac=json.loads((F/'converted'/stem/'ACCEPTANCE.json').read_text());assert ac['acquisition_pass'];assert (F/'raw'/name).stat().st_size==item['bytes'];update(name,stage='reused_complete');continue
  raw=F/'raw'/name;part=F/'raw'/(name+'.server.partial');upload=F/'raw'/(name+'.uploading');assert not raw.exists() and not upload.exists(),'Existing incomplete processing requires audit'
  attempts=[];ok=False;stop_source=False
  for attempt in range(2):
   before=part.stat().st_size if part.exists() else 0;assert before<=item['bytes'];remaining=item['bytes']-before;seconds=max(1800,min(43200,math.ceil(remaining/8192*1.5)+300));update(name,stage='downloading',attempt=attempt+1,prefix_bytes=before,expected_bytes=item['bytes'],deadline_seconds=seconds)
   cmd=['curl','-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time',str(seconds),'--speed-limit','1024','--speed-time','180','--silent','--show-error']
   if before:cmd+=['--continue-at','-']
   cmd+=['--output',str(part),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name]
   r=subprocess.run(cmd,capture_output=True,text=True);after=part.stat().st_size if part.exists() else 0;attempts.append({'exit':r.returncode,'before':before,'after':after,'stderr':r.stderr,'utc':now()});save(W/(name+'.transfer.json'),attempts)
   if r.returncode==0 and after==item['bytes']:ok=True;break
   if r.returncode in [5,6,7,60]:stop_source=True;break
   if r.returncode not in [28,56] or after<=before:break
  if not ok:
   update(name,stage='transfer_incomplete',attempts=attempts)
   if stop_source:state['stop_reason']='source connectivity or certificate failure; no repeated queue attempts';break
   continue
  if name in prefixes:
   size,digest=prefixes[name];assert sha(part,size)==digest
  with part.open('rb') as h:assert h.read(2)==b'\x01\xa1'
  digest=sha(part);save(W/(name+'.raw_receipt.json'),{'bytes':item['bytes'],'sha256':digest,'transfer_exit':0,'certificate_validation_enabled':True});assert not raw.exists() and not upload.exists();part.rename(upload)
  update(name,stage='conversion_and_search',sha256=digest)
  with (W/(name+'.pipeline.log')).open('w',encoding='utf-8') as log:r=subprocess.run(['python3',str(P/'scripts/cgnl_c04_accept_convert_file.py'),name,digest],stdout=log,stderr=subprocess.STDOUT)
  if r.returncode!=0 or not complete.exists():update(name,stage='pipeline_failed',exit=r.returncode);continue
  update(name,stage='accepted',raw_sha256=digest,complete_receipt=str(complete))
 state['complete']=len(state['files'])==13 and all(x['stage'] in ['accepted','reused_complete','reused_pilot'] for x in state['files'].values())
 if state['complete']:
  r=subprocess.run(['python3',str(P/'scripts/cgnl_c04_joint_psm_q.py')],capture_output=True,text=True);z={'joint_exit':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
  if r.returncode==0:
   q=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(P/'scripts/cgnl_c04_joint_q_independent.R'),str(F/'joint_identification')],capture_output=True,text=True);z.update(R_exit=q.returncode,R_stdout=q.stdout,R_stderr=q.stderr)
  save(W/'JOINT_RUN_EXIT.json',z)
except Exception as e:state['error']=repr(e)
finally:
 state['finished']=now();save(W/'STATUS.json',state);save(W/'EXIT.json',{'acquisition_search_complete':state['complete'],'error':state.get('error'),'utc':now(),'whole_project_complete':False});print(json.dumps(state),flush=True)
raise SystemExit(0 if state['complete'] else 1)
