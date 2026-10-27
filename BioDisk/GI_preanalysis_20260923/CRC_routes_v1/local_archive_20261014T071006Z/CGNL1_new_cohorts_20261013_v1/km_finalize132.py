from pathlib import Path
import subprocess,json,os,re,gzip,hashlib,datetime
R=Path('/data/nas1/wangmingyuan_OD/CGNL1_public_download_20261013_v1')
ps=subprocess.run(['ps','-u',str(os.getuid()),'-o','pid=,args='],capture_output=True,text=True,check=True).stdout
for line in ps.splitlines():
 pid,cmd=line.strip().split(None,1)
 if 'python' in cmd and cmd.endswith('/CGNL1_public_download_20261013_v1/km_segment132.py'):subprocess.run(['tail','--pid='+pid,'-f','/dev/null'])
name='GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz';dest=R/name
if dest.exists() and (R/'SEGMENT132_COMPLETE.json').exists():print('ALREADY_VERIFIED');raise SystemExit(0)
plan=json.loads((R/'SEGMENT132_PLAN.json').read_text());D=R/'GSE132465_segments';url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/'+name
attempts=[]
for i,start,end in plan['chunks']:
 p=D/(str(i)+'.part');expected=end-start+1
 for a in range(3):
  have=p.stat().st_size if p.exists() else 0
  if have==expected:break
  assert have<expected
  add=D/(str(i)+'.retry'+str(a));hdr=D/(str(i)+'.retry'+str(a)+'.headers');log=D/(str(i)+'.retry'+str(a)+'.log')
  with log.open('w') as f:
   proc=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','300','--range',str(start+have)+'-'+str(end),'--dump-header',str(hdr),'--output',str(add),url],stderr=f)
  header=hdr.read_text().lower();ok=('content-range: bytes '+str(start+have)+'-'+str(end)+'/'+str(plan['total'])) in header
  attempts.append(dict(segment=i,attempt=a,exit_code=proc.returncode,range_start=start+have,expected_end=end,bytes=add.stat().st_size if add.exists() else 0,header_valid=ok))
  if ok and add.exists() and add.stat().st_size<=expected-have:
   with p.open('ab') as f:f.write(add.read_bytes())
 assert p.stat().st_size==expected,'Segment remains incomplete '+str(i)
(R/'SEGMENT132_RECOVERY.json').write_text(json.dumps(attempts,indent=2))
prefix=R/(name+'.part');assert prefix.stat().st_size==plan['prefix_bytes']
assembled=R/(name+'.assembled')
with assembled.open('wb') as f:
 for p in [prefix]+[D/(str(x[0])+'.part') for x in plan['chunks']]:
  with p.open('rb') as inp:
   for b in iter(lambda:inp.read(1048576),b''):f.write(b)
assert assembled.stat().st_size==plan['total']
with gzip.open(assembled,'rb') as f:
 while f.read(1048576):pass
assembled.replace(dest)
rec=dict(dataset='GSE132465',file=name,url=url,status='verified',bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'SEGMENT132_COMPLETE.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
