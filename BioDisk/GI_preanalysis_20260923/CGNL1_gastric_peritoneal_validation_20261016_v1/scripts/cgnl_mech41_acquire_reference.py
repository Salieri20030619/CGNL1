from pathlib import Path
import requests,json,hashlib,gzip,datetime,concurrent.futures,time
O=Path.cwd()/'mechanism_investigation_v2';R=O/'reference/gencode_v47';R.mkdir(parents=True,exist_ok=True);base='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/'
names=['MD5SUMS','gencode.v47.transcripts.fa.gz','GRCh38.primary_assembly.genome.fa.gz']
def get(n):
 p=R/n;tmp=R/(n+'.partial');r={'name':n,'url':base+n,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  assert not p.exists() and not tmp.exists(),'existing input; inspect before restart'
  md5=hashlib.md5();sha=hashlib.sha256();size=0;start=time.monotonic()
  with requests.get(base+n,stream=True,timeout=(20,120)) as a:
   a.raise_for_status();r['content_length']=a.headers.get('Content-Length')
   with tmp.open('wb') as f:
    for b in a.iter_content(1024*1024):
     if b:f.write(b);md5.update(b);sha.update(b);size+=len(b)
     if time.monotonic()-start>14400:raise TimeoutError('Four-hour transfer cap')
  if r['content_length']:assert size==int(r['content_length'])
  if n.endswith('.gz'):
   with gzip.open(tmp,'rb') as f:
    assert f.read(1)==b'>'
    while f.read(4*1024*1024):pass
  tmp.rename(p);r.update(success=True,bytes=size,md5=md5.hexdigest(),sha256=sha.hexdigest())
 except Exception as e:r.update(success=False,error=repr(e))
 r['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'metadata'/('REFERENCE_'+n.replace('.','_')+'.json')).write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:rs=list(ex.map(get,names))
(O/'metadata/REFERENCE_DOWNLOAD_COMPLETE.json').write_text(json.dumps(rs,indent=2),encoding='utf-8')
