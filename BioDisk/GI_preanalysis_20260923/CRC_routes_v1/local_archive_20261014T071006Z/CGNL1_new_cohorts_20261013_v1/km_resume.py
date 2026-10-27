from pathlib import Path
import subprocess,os,json,hashlib,gzip,datetime,concurrent.futures
R=Path('/data/nas1/wangmingyuan_OD/CGNL1_public_download_20261013_v1')
ps=subprocess.run(['ps','-u',str(os.getuid()),'-o','pid=,args='],capture_output=True,text=True,check=True).stdout
for line in ps.splitlines():
 pid,cmd=line.strip().split(None,1)
 if 'python' in cmd and cmd.endswith('/CGNL1_public_download_20261013_v1/km_download.py'):
  subprocess.run(['tail','--pid='+pid,'-f','/dev/null'],check=True)
if (R/'MANIFEST.json').exists():(R/'MANIFEST_first_attempt.json').write_bytes((R/'MANIFEST.json').read_bytes())
def fetch(item):
 ds,ext=item;name=ds+'_GEO_processed_CRC_10X_raw_UMI_count_matrix.'+ext+'.gz'
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/'+ds[:-3]+'nnn/'+ds+'/suppl/'+name
 dest=R/name;tmp=R/(name+'.part');log=R/(ds+'_resume.log')
 rc=0
 if not dest.exists():
  with log.open('w') as f:
   p=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','900','--continue-at','-','--speed-time','90','--speed-limit','2048','--output',str(tmp),url],stdout=f,stderr=subprocess.STDOUT)
  rc=p.returncode
  if rc==0:
   with gzip.open(tmp,'rb') as f:
    while f.read(1048576):pass
   tmp.replace(dest)
 rec=dict(dataset=ds,file=name,url=url,exit_code=rc,status='verified' if dest.exists() else 'failed')
 if dest.exists():rec.update(bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
 else:rec.update(log=log.read_text())
 print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:out=list(p.map(fetch,[('GSE200997','csv'),('GSE132465','txt')]))
(R/'MANIFEST.json').write_text(json.dumps(out,indent=2));(R/'COMPLETE.json').write_text(json.dumps(dict(all_success=all(x['status']=='verified' for x in out),utc=datetime.datetime.now(datetime.timezone.utc).isoformat())))
