import os,json,time,hashlib,gzip,urllib.request,urllib.parse,re,shutil,traceback
from pathlib import Path
BASE=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
ROOT=BASE/'data_prepare_20260929'
ROOT.mkdir(exist_ok=True)
OP=urllib.request.build_opener(urllib.request.ProxyHandler({}))
RESERVE=20*1024**3

def atomic(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def fetch(u):
 with OP.open(urllib.request.Request(u,headers={'User-Agent':'CGNL1-research-data-preparation/1.0'}),timeout=45) as r:return r.read()
def verify(p):
 if p.stat().st_size==0:raise ValueError('empty file')
 with p.open('rb') as f:head=f.read(128)
 if head.lstrip().lower().startswith((b'<!doctype html',b'<html')):raise ValueError('HTML returned instead of data')
 if p.name.endswith('.gz'):
  with gzip.open(p,'rb') as f:
   while f.read(4*1024**2):pass
 if p.name.endswith(('.h5seurat','.h5ad','.h5')) and not head.startswith(b'\x89HDF\r\n\x1a\n'):raise ValueError('not HDF5')
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4*1024**2),b''):h.update(block)
 return h.hexdigest()

def run(wave,jobs):
 dest=ROOT/wave;dest.mkdir(exist_ok=True);atomic(dest/'PLAN.json',jobs)
 results=[]
 for i,j in enumerate(jobs):
  p=dest/j['dataset']/j['name'];p.parent.mkdir(exist_ok=True)
  rec=dict(j,path=str(p));start=time.time()
  try:
   if p.exists():
    rec.update(status='verified_existing',bytes=p.stat().st_size,sha256=verify(p));results.append(rec);continue
   part=Path(str(p)+'.part')
   for attempt in range(1,5):
    try:
     offset=part.stat().st_size if part.exists() else 0
     headers={'User-Agent':'CGNL1-research-data-preparation/1.0'}
     if offset:headers['Range']=f'bytes={offset}-'
     with OP.open(urllib.request.Request(j['url'],headers=headers),timeout=60) as r:
      if offset and r.status!=206:offset=0
      remaining=int(r.headers.get('Content-Length','0'))
      if shutil.disk_usage(ROOT).free-max(remaining,1024**3)<RESERVE:raise RuntimeError('disk_reserve_20GiB')
      total=offset+remaining if remaining else None
      with part.open('ab' if offset else 'wb') as f:
       n=offset;last=0
       while True:
        chunk=r.read(1024**2)
        if not chunk:break
        f.write(chunk);n+=len(chunk)
        if time.time()-last>20:
         if shutil.disk_usage(ROOT).free<RESERVE:raise RuntimeError('disk_reserve_20GiB')
         atomic(dest/'PROGRESS.json',{'index':i+1,'jobs':len(jobs),'name':j['name'],'bytes':n,'expected':total,'time':time.time()});last=time.time()
      if total is not None and n!=total:raise ValueError(f'length mismatch {n}/{total}')
     # Verify using the final suffix before marking complete.
     part.replace(p)
     try:sha=verify(p)
     except Exception:
      p.replace(part);raise
     if j.get('sha256') and sha.lower()!=j['sha256'].lower():
      p.replace(part);raise ValueError('publisher sha256 mismatch')
     rec.update(status='downloaded_verified',bytes=p.stat().st_size,sha256=sha,attempt=attempt);break
    except Exception as e:
     rec['last_error']=str(e)
     if 'disk_reserve' in str(e) or attempt==4:raise
     time.sleep(min(10*attempt,30))
  except Exception as e:rec.update(status='blocked' if 'disk_reserve' in str(e) else 'failed',error=str(e))
  rec['seconds']=round(time.time()-start,2);results.append(rec)
  atomic(dest/'MANIFEST.json',results);print(json.dumps(rec,ensure_ascii=False),flush=True)
  if rec['status']=='blocked':break
 atomic(dest/'MANIFEST.json',results)
 code=0 if len(results)==len(jobs) and all(x['status'] in ['downloaded_verified','verified_existing'] for x in results) else 1
 atomic(dest/'COMPLETE.json',{'exit_code':code,'planned':len(jobs),'recorded':len(results),'verified':sum(x['status'] in ['downloaded_verified','verified_existing'] for x in results),'time':time.time()})
 return code

if __name__=='__main__':
 import sys
 if len(sys.argv)>1:
  plan=json.loads(Path(sys.argv[1]).read_text());sys.exit(run(Path(sys.argv[1]).stem,plan))
 wave='wave01_core';jobs=[]
 m=BASE/'expanded_screen_20260925/metadata'
 soft=gzip.decompress((m/'GSE234804.soft.gz').read_bytes()).decode()
 for line in soft.splitlines():
  if line.startswith('!Sample_supplementary_file'):
   u=line.split(' = ',1)[1].replace('ftp://','https://')
   if u!='NONE':jobs.append({'dataset':'GSE234804','name':u.rsplit('/',1)[1],'url':u,'role':'human_scRNA_processed'})
 for g in ['GSE131418']:
  u=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/matrix/'
  raw=fetch(u);(ROOT/f'{g}_matrix_index.html').write_bytes(raw)
  links=re.findall(r'href="([^"]+series_matrix[^"?]+\.gz)"',raw.decode())
  for name in sorted(set(links)):
   jobs.append({'dataset':g,'name':name.rsplit('/',1)[1],'url':urllib.parse.urljoin(u,name),'role':'bulk_expression_with_metadata'})
 atomic(ROOT/'wave01_core.json',jobs)
 sys.exit(run(wave,jobs))
