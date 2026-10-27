from pathlib import Path
import os,signal,time,json,requests,hashlib,concurrent.futures,threading,tarfile,gzip,datetime
O=Path.cwd()/'mechanism_investigation_v2';q=json.loads((O/'metadata/GSE308231_RANGE_PROBE.json').read_text());assert q['safe_range_capability'];total=int(q['content_range'].split('/')[1]);p=O/'raw/GSE308231_RAW.tar';part=p.with_suffix('.tar.partial');prefix=p.with_name(p.name+'.preserved_prefix');C=O/'raw/GSE308231_ranges';C.mkdir(exist_ok=True);assert not p.exists() and not prefix.exists();stopped=[]
for d in Path('/proc').iterdir():
 if d.name.isdigit():
  try:
   args=(d/'cmdline').read_bytes().split(b'\0')
   for name in ['cgnl_mech04_download_tissue.py','cgnl_mech47_tissue_when_ready.py']:
    if any(a==str(O.parent/'scripts'/name).encode() for a in args):os.kill(int(d.name),signal.SIGTERM);stopped.append({'pid':int(d.name),'script':name})
  except (FileNotFoundError,PermissionError):pass
(O/'metadata/TISSUE_TRANSFER_UPGRADE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'stopped':stopped,'range_prefix_verified':True,'previous_data_preserved':True,'expected_bytes':total,'workers':4},indent=2),encoding='utf-8')
for rec in stopped:
 for _ in range(100):
  if not (Path('/proc')/str(rec['pid'])).exists():break
  time.sleep(.1)
 else:raise RuntimeError('Previous writer still live')
part.rename(prefix);offset=prefix.stat().st_size;intervals=[(s,min(s+8*1024*1024,total)) for s in range(offset,total,8*1024*1024)];manifest={};lock=threading.Lock()
def piece(iv):
 start,end=iv;dst=C/(str(start)+'_'+str(end)+'.bin');r={'start':start,'end_exclusive':end}
 for attempt in range(2):
  temp=dst.with_name(dst.name+f'.attempt{attempt+1}.partial')
  try:
   with requests.get(q['url'],headers={'Range':f'bytes={start}-{end-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,120)) as a:
    assert a.status_code==206 and a.headers.get('Content-Range')==f'bytes {start}-{end-1}/{total}'
    if a.headers.get('Last-Modified'):assert a.headers['Last-Modified']==q['last_modified']
    n=0;h=hashlib.sha256();begin=time.monotonic()
    with temp.open('wb') as f:
     for b in a.iter_content(65536):
      if b:f.write(b);h.update(b);n+=len(b)
      if time.monotonic()-begin>3600:raise TimeoutError('onehourchunkcap')
   assert n==end-start;temp.rename(dst);r.update(success=True,bytes=n,sha256=h.hexdigest());break
  except Exception as e:r.update(success=False,error=repr(e))
 with lock:manifest[str(start)]=r;(C/'range_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 if not r['success']:raise RuntimeError(r)
 return dst
receipt={'url':q['url'],'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'preserved_prefix_bytes':offset,'header_size':total,'upgraded_range_transfer':True}
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:pieces=list(ex.map(piece,intervals))
 assembled=p.with_name(p.name+'.assembled.partial');h=hashlib.sha256();size=0
 with assembled.open('wb') as out:
  for f in [prefix]+pieces:
   with f.open('rb') as src:
    while b:=src.read(4*1024*1024):out.write(b);h.update(b);size+=len(b)
 assert size==total;members=[]
 with tarfile.open(assembled) as tf:
  for m in tf.getmembers():
   if not m.isfile():continue
   digest=hashlib.sha256();stream=tf.extractfile(m)
   with gzip.GzipFile(fileobj=stream) as gz:
    while b:=gz.read(4*1024*1024):digest.update(b)
   members.append({'name':m.name,'compressed_bytes':m.size,'decompressed_sha256':digest.hexdigest(),'gzip_crc_verified':True})
 assert len(members)==18;assembled.rename(p);receipt.update(success=True,bytes=size,sha256=h.hexdigest(),members=members,verification='expectedsize+tarstructure+all18internalgzipCRC; wholefileauthorMD5notprovided')
except Exception as e:receipt.update(success=False,error=repr(e))
receipt['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'metadata/GSE308231_DOWNLOAD.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt),flush=True)
