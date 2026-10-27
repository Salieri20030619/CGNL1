from pathlib import Path
import subprocess,concurrent.futures,json,gzip,hashlib,datetime
R=Path('/data/nas1/wangmingyuan_OD/CGNL1_public_download_20261013_v1');R.mkdir(exist_ok=True)
def fetch(item):
 ds,ext=item;name=ds+'_GEO_processed_CRC_10X_raw_UMI_count_matrix.'+ext+'.gz'
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/'+ds[:-3]+'nnn/'+ds+'/suppl/'+name
 dest=R/name;tmp=R/(name+'.part');log=R/(ds+'.log')
 with log.open('w') as f:
  p=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','300','--speed-time','45','--speed-limit','2048','--output',str(tmp),url],stdout=f,stderr=subprocess.STDOUT)
 rec=dict(dataset=ds,file=name,url=url,exit_code=p.returncode)
 if p.returncode==0:
  with gzip.open(tmp,'rb') as f:
   while f.read(1048576):pass
  tmp.replace(dest);rec.update(status='verified',bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
 else:rec.update(status='failed',log=log.read_text())
 print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:out=list(p.map(fetch,[('GSE200997','csv'),('GSE132465','txt')]))
(R/'MANIFEST.json').write_text(json.dumps(out,indent=2))
(R/'COMPLETE.json').write_text(json.dumps(dict(all_success=all(x['status']=='verified' for x in out),utc=datetime.datetime.now(datetime.timezone.utc).isoformat())))
