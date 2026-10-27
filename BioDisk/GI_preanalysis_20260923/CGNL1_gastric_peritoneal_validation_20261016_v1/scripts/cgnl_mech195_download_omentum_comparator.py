
from pathlib import Path
import requests,json,datetime,gzip,hashlib,time,concurrent.futures,threading
p=Path.cwd();m=p/'mechanism_investigation_v2';dest=m/'raw/GSE286057_comparator_195';dest.mkdir(exist_ok=True)
spec=json.loads((m/'GSE286057_COMPARATOR_SPECIFICATION_195.json').read_text(encoding='utf-8'));tasks=[(s['accession'],url) for s in spec['samples'] for url in s['urls']];tls=threading.local()
def one(task):
 acc,url=task;name=url.rsplit('/',1)[1];f=dest/name;rec={'sample':acc,'url':url,'file':str(f)}
 if not hasattr(tls,'s'):tls.s=requests.Session();tls.s.trust_env=False
 for attempt in range(1,4):
  tmp=f.with_name(f.name+'.attempt'+str(attempt)+'.partial')
  try:
   assert not f.exists(),'Destination already exists, review before repeat'
   with tls.s.get(url,stream=True,timeout=(20,90)) as r:
    r.raise_for_status();expected=int(r.headers['Content-Length']) if 'Content-Length' in r.headers else None
    h=hashlib.sha256();n=0;start=time.monotonic()
    with tmp.open('xb') as out:
     for b in r.iter_content(262144):
      out.write(b);h.update(b);n+=len(b)
      if time.monotonic()-start>7200:raise TimeoutError('bounded two-hour file transfer')
   assert expected is None or expected==n
   with gzip.open(tmp,'rb') as inp:
    while inp.read(1048576):pass
   tmp.rename(f);rec.update(success=True,bytes=n,content_length=expected,sha256=h.hexdigest(),gzip_crc_verified=True,attempt=attempt);return rec
  except Exception as e:
   rec.update(success=False,error=repr(e),attempt=attempt)
   if attempt<3:time.sleep(30*attempt)
 return rec
receipts=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for r in ex.map(one,tasks):
  receipts.append(r)
  (m/'metadata/GSE286057_DOWNLOAD_RECEIPTS_195.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
  print(r['sample'],Path(r['file']).name,r['success'],flush=True)
complete={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':len(receipts)==108 and all(r['success'] for r in receipts),'files':len(receipts),'source':'official GEO per-sample processed matrices; no raw human reads accessed'}
(m/'metadata/GSE286057_DOWNLOAD_COMPLETE_195.json').write_text(json.dumps(complete,indent=2),encoding='utf-8')
assert complete['all_pass'],complete
