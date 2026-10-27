from pathlib import Path
import subprocess,json,time,datetime,hashlib,zipfile,concurrent.futures,os,traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=B/'raw/intermittent_564';D.mkdir(exist_ok=True)
M=B/'metadata/intermittent_564';M.mkdir(exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,obj):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2),encoding='utf-8');tmp.replace(path)
jobs=[dict(id='HCC',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10783603/supplementaryFiles',kind='zip',wanted=['mmc2.xlsx','mmc3.xlsx','mmc6.xlsx']),dict(id='CRC',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12686132/supplementaryFiles',kind='zip',wanted=['MCO2-6-e70492-s011.xlsx']),dict(id='ESCC',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8473482/supplementaryFiles',kind='zip',wanted=['CTM2-11-e538-s001.doc']),dict(id='TAO',url='https://www.rsc.org/suppdata/d0/mo/d0mo00125b/d0mo00125b1.zip',kind='zip',wanted=None),dict(id='CCA',url='https://zenodo.org/api/records/7708351/files/Supplementary_Table_1.xlsx/content',kind='xlsx',md5='0e36788a46cfb335dc9562844d953528',size=68506395)]
policy=dict(single_transfer_timeout_seconds=3600,connect_timeout_seconds=60,retry_interval_seconds=900,max_attempts=3,concurrency=2,low_speed_limit_bytes_per_second=1024,low_speed_time_seconds=300,max_response_bytes=262144000,transport='direct; no active proxy listener observed; do not change system proxy',resume='Static CCA file supports validated full-file MD5; generated supplement ZIPs use separate fresh attempts because byte stability is not established',access_denial='HTTP401/403 retry once after900s then stop; no credential or access-control changes')
write(M/'PLAN.json',dict(created=now(),pid=os.getpid(),policy=policy,jobs=jobs))
def run(job):
 state=dict(id=job['id'],url=job['url'],stage='queued',attempts=[],utc=now());path=M/(job['id']+'.json');write(path,state)
 try:
  for attempt in range(1,4):
   dest=D/(job['id']+('_static.part' if job['kind']=='xlsx' else '_attempt'+str(attempt)+'.zip.part'))
   if job['id']=='CCA' and attempt==1 and not dest.exists():
    prior=B/'raw/qualification_558/CCA_Table1.xlsx'
    if prior.exists() and prior.stat().st_size<job['size']:
     import shutil
     shutil.copyfile(prior,dest)
   state.update(stage='downloading',utc=now(),attempt=attempt,current_file=str(dest));write(path,state)
   cmd=['curl','-q','--noproxy','*','--location','--fail','--silent','--show-error','--connect-timeout','60','--max-time','3600','--speed-limit','1024','--speed-time','300','--max-filesize','262144000','--write-out','%{http_code}','--output',str(dest)]
   if job['kind']=='xlsx' and dest.exists() and dest.stat().st_size:cmd+=['--continue-at','-']
   r=subprocess.run(cmd+[job['url']],capture_output=True,text=True);status=r.stdout.strip();entry=dict(attempt=attempt,utc=now(),exit=r.returncode,http=status,error=r.stderr[-3000:],bytes=dest.stat().st_size if dest.exists() else 0);state['attempts'].append(entry)
   if r.returncode==0:
    state.update(stage='validating',utc=now());write(path,state)
    if job.get('size') is not None:assert dest.stat().st_size==job['size'],'Size mismatch'
    if job.get('md5'):
     md=hashlib.md5()
     with dest.open('rb') as f:
      for chunk in iter(lambda:f.read(8*1024*1024),b''):md.update(chunk)
     assert md.hexdigest()==job['md5'],'Publisher MD5 mismatch';state['md5']=md.hexdigest()
    with zipfile.ZipFile(dest) as z:
     assert z.testzip() is None,'ZIP CRC failure';names=z.namelist();wanted=job.get('wanted');chosen=[n for n in names if wanted is None or Path(n).name in wanted]
     if wanted is not None:assert set(wanted)<=set(Path(n).name for n in chosen),'Required members missing'
     extracted=[]
     if job['kind']=='zip':
      edir=D/(job['id']+'_verified');edir.mkdir(exist_ok=True)
      for n in chosen:
       if n.endswith('/'):continue
       if wanted is None and not n.lower().endswith(('.xlsx','.xls','.csv','.tsv','.doc','.docx','.pdf','.txt')):continue
       target=edir/Path(n).name
       assert not target.exists(),'Refuse overwrite of verified member'
       with z.open(n) as source,target.open('wb') as target_f:
        import shutil
        shutil.copyfileobj(source,target_f)
       extracted.append(str(target))
    final=dest.with_suffix('.verified');assert not final.exists();dest.rename(final);state.update(stage='complete_verified',utc=now(),file=str(final),zip_crc_pass=True,extracted=extracted,scientific_qualification='pending; download completion is not cohort validation');write(path,state);return state
   if r.returncode in [33,63] or status=='404':state['terminal_reason']='Range/size/not-found requires manual source review';break
   if status in ['401','403'] and attempt>=2:state['terminal_reason']='Persistent access denial after one delayed retry';break
   if r.returncode not in [7,18,22,28,35,52,55,56]:state['terminal_reason']='Non-transient error';break
   if attempt<3:
    state.update(stage='waiting_retry',utc=now(),next_retry_utc=datetime.datetime.fromtimestamp(time.time()+900,datetime.timezone.utc).isoformat());write(path,state);time.sleep(900)
  state.update(stage='stopped_unverified',utc=now());write(path,state);return state
 except Exception as e:
  state.update(stage='stopped_unverified',utc=now(),error=repr(e),traceback=traceback.format_exc());write(path,state);return state
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
write(M/'EXIT.json',dict(utc=now(),all_verified=all(r['stage']=='complete_verified' for r in results),results=[dict(id=r['id'],stage=r['stage']) for r in results]))
