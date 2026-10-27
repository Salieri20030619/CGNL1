from pathlib import Path
import json,hashlib,requests,time,concurrent.futures,datetime,os
b=Path(__file__).resolve().parents[1];m=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text());dest=Path(m['destination']);chunk_size=4194304
def digest(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for data in iter(lambda:f.read(8388608),b''):h.update(data)
 return h.hexdigest()
def fetch(item,start,end):
 for attempt in range(5):
  try:
   with requests.Session() as session:
    session.trust_env=False
    with session.get(item['url'],headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(15,90)) as r:
     r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{item["bytes"]}'
     data=r.raw.read(end-start+1);assert len(data)==end-start+1;return data
  except Exception:
   if attempt==4:raise
   time.sleep(5)
def run(item):
 target=dest/item['name'];part=dest/(item['name']+'.part');total=item['bytes']
 if target.exists():assert target.stat().st_size==total and digest(target)==item['md5'];return dict(item,verified=True,local_path=str(target))
 offset=part.stat().st_size if part.exists() else 0;assert offset<=total
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool,part.open('ab') as f:
  assert f.tell()==offset;pending=[];next_offset=offset
  def enqueue():
   nonlocal next_offset
   if next_offset<total:
    end=min(next_offset+chunk_size,total)-1;pending.append((next_offset,end,pool.submit(fetch,item,next_offset,end)));next_offset=end+1
  enqueue();enqueue()
  while pending:
   start,end,future=pending.pop(0);data=future.result();assert f.tell()==start and len(data)==end-start+1
   f.write(data);f.flush();os.fsync(f.fileno());assert f.tell()==end+1
   p=b/'metadata'/('DOWNLOAD468_'+item['name']+'.json');tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':end+1,'total':total,'pid':os.getpid(),'server_direct':True,'connections_per_file':2}),encoding='utf-8');tmp.replace(p);enqueue()
 assert part.stat().st_size==total;h=digest(part);assert h==item['md5'],'MD5 mismatch; part retained'
 part.rename(target);return dict(item,verified=True,local_path=str(target))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:receipts=list(pool.map(run,m['items']))
(b/'metadata/PILOT_DOWNLOAD_COMPLETE_468.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'run':m['run'],'files':receipts,'all_md5_verified':True,'server_direct':True},indent=2),encoding='utf-8');print('COMPLETE',flush=True)
