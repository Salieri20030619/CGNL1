from pathlib import Path
import concurrent.futures, subprocess, json, time, os, hashlib, zipfile, traceback, fcntl
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H02')
M=B/'metadata/resume_581'; M.mkdir(exist_ok=True)
R=D/'resume_581'; R.mkdir(exist_ok=True)
URL='https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD008373/Phospho_MaxQuant_results.zip'
N=17879997918
lock=(R/'lock').open('w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def save(name,obj):
 p=M/name; q=p.with_suffix('.tmp'); q.write_text(json.dumps(obj,indent=2),encoding='utf-8');q.replace(p)
def fetch(a,b,p):
 for attempt in range(2):
  h=p.with_suffix('.headers')
  r=subprocess.run(['curl','--proxy','http://127.0.0.1:7897','--noproxy','','--fail','--silent','--show-error','--connect-timeout','20','--max-time','900','--range',f'{a}-{b}','-D',str(h),'-o',str(p),URL],capture_output=True,text=True)
  headers=h.read_text().lower() if h.exists() else ''
  if r.returncode==0 and p.stat().st_size==b-a+1 and f'content-range: bytes {a}-{b}/{N}' in headers:return
  if attempt==1:raise RuntimeError(str((a,b,r.returncode,r.stderr[-500:],headers[-1000:])))
  time.sleep(5)
started=time.time();save('RUN.json',dict(pid=os.getpid(),started=started,url=URL,expected_bytes=N))
try:
 final=D/'Phospho_MaxQuant_results.zip'
 if final.exists():raise RuntimeError('Existing final archive; refusing overwrite')
 seed=max(D.glob('Phospho_MaxQuant_results.zip.attempt*.part'),key=lambda p:p.stat().st_size)
 size=seed.stat().st_size;assert 0<size<N
 with seed.open('rb') as f:
  for a,b in [(0,1048575),(size-1048576,size-1)]:
   p=R/f'probe_{a}.bin';fetch(a,b,p);f.seek(a);assert f.read(b-a+1)==p.read_bytes(),'Existing prefix mismatch'
 save('PREFIX_QA.json',dict(seed=str(seed),bytes=size,head_and_tail_source_match=True))
 step=128*1024*1024
 chunks=[(a,min(a+step,N)-1,R/f'{a:012d}.part') for a in range(size,N,step)]
 done=0
 def worker(t):
  a,b,p=t; receipt=p.with_suffix('.verified.json')
  if receipt.exists() and p.exists() and p.stat().st_size==b-a+1:
   old=json.loads(receipt.read_text())
   if hashlib.sha256(p.read_bytes()).hexdigest()==old['sha256']:return t
  fetch(a,b,p)
  receipt.write_text(json.dumps(dict(start=a,end=b,sha256=hashlib.sha256(p.read_bytes()).hexdigest())),encoding='utf-8')
  return t
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  for a,b,p in pool.map(worker,chunks):
   done+=b-a+1;save('PROGRESS.json',dict(stage='downloading',seed_bytes=size,downloaded_new_bytes=done,total=N,elapsed=time.time()-started))
 part=R/'assembled.zip.partial';sha=hashlib.sha256()
 with part.open('wb') as out:
  for p in [seed]+[t[2] for t in chunks]:
   with p.open('rb') as f:
    for buf in iter(lambda:f.read(8*1024*1024),b''):out.write(buf);sha.update(buf)
 assert part.stat().st_size==N
 save('PROGRESS.json',dict(stage='zip_crc_verification',bytes=N,sha256=sha.hexdigest()))
 with zipfile.ZipFile(part) as z:
  assert z.testzip() is None,'ZIP CRC failed'
  members=[dict(name=i.filename,size=i.file_size) for i in z.infolist()]
  dest=D/'selected_581';dest.mkdir(exist_ok=True);selected=[]
  for i in z.infolist():
   if Path(i.filename).name.lower() not in ['phospho (sty)sites.txt','proteingroups.txt','summary.txt','parameters.txt']:continue
   target=(dest/i.filename).resolve();assert dest.resolve() in target.parents
   assert (i.external_attr>>16)&0o170000!=0o120000
   z.extract(i,dest);selected.append(str(target))
 part.rename(final)
 receipt=dict(stage='verified_complete',bytes=N,sha256=sha.hexdigest(),zip_crc_pass=True,members=members,selected=selected,url=URL,elapsed=time.time()-started)
 save('DOWNLOAD_COMPLETE.json',receipt)
 final.with_suffix('.zip.receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 save('EXIT.json',dict(exit_code=0,stage='verified_complete',utc=time.time()))
 print(json.dumps(receipt),flush=True)
except Exception:
 save('EXIT.json',dict(exit_code=1,traceback=traceback.format_exc(),utc=time.time()))
 raise
