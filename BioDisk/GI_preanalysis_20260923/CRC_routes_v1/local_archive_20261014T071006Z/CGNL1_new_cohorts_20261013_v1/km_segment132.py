from pathlib import Path
import subprocess,json,os,signal,concurrent.futures,re,gzip,hashlib,datetime
R=Path('/data/nas1/wangmingyuan_OD/CGNL1_public_download_20261013_v1')
name='GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz';dest=R/name;prefix=R/(name+'.part')
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/'+name
# Stop only the curl child writing this task's exact prefix file.
ps=subprocess.run(['ps','-u',str(os.getuid()),'-o','pid=,args='],capture_output=True,text=True,check=True).stdout
stopped=[]
for line in ps.splitlines():
 pid,cmd=line.strip().split(None,1)
 if cmd.startswith('curl ') and str(prefix) in cmd and url in cmd:
  os.kill(int(pid),signal.SIGTERM);subprocess.run(['tail','--pid='+pid,'-f','/dev/null']);stopped.append(int(pid))
if dest.exists():print('ALREADY_COMPLETE');raise SystemExit(0)
D=R/'GSE132465_segments';D.mkdir(exist_ok=True)
head=D/'probe.headers';probe=D/'probe.byte'
p=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','60','--range','0-0','--dump-header',str(head),'--output',str(probe),url])
assert p.returncode==0
h=head.read_text().lower();total=int(re.findall(r'content-range:\s*bytes 0-0/(\d+)',h)[-1]);start=prefix.stat().st_size
assert 0<start<total
step=(total-start+2)//3
chunks=[(i,s,min(s+step,total)-1) for i,s in enumerate(range(start,total,step))]
(R/'SEGMENT132_PLAN.json').write_text(json.dumps(dict(stopped_curl=stopped,prefix_bytes=start,total=total,chunks=chunks),indent=2))
def get(rec):
 i,s,e=rec;out=D/(str(i)+'.part');hdr=D/(str(i)+'.headers');log=D/(str(i)+'.log')
 with log.open('w') as err:
  p=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','900','--speed-time','90','--speed-limit','2048','--range',str(s)+'-'+str(e),'--dump-header',str(hdr),'--output',str(out),url],stderr=err)
 ok=p.returncode==0 and out.stat().st_size==e-s+1 and ('content-range: bytes '+str(s)+'-'+str(e)+'/'+str(total)) in hdr.read_text().lower()
 print('SEGMENT',i,'EXIT',p.returncode,'OK',ok,flush=True)
 return dict(part=i,path=str(out),exit_code=p.returncode,verified=ok,start=s,end=e)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:parts=list(pool.map(get,chunks))
(R/'SEGMENT132_PARTS.json').write_text(json.dumps(parts,indent=2));assert all(p['verified'] for p in parts)
assembled=R/(name+'.assembled')
with assembled.open('wb') as f:
 for source in [prefix]+[Path(p['path']) for p in parts]:
  with source.open('rb') as inp:
   for b in iter(lambda:inp.read(1048576),b''):f.write(b)
assert assembled.stat().st_size==total
with gzip.open(assembled,'rb') as f:
 while f.read(1048576):pass
assembled.replace(dest)
rec=dict(dataset='GSE132465',file=name,url=url,status='verified',bytes=total,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'SEGMENT132_COMPLETE.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
