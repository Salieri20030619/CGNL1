from pathlib import Path
import requests,json,hashlib,time,datetime,concurrent.futures,os
b=Path(__file__).resolve().parents[1];m=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text());dest=Path(m['destination'])
def digest(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8388608),b''):h.update(chunk)
 return h.hexdigest()
def download(item):
 session=requests.Session();session.trust_env=False
 target=dest/item['name'];part=dest/(item['name']+'.part');progress=b/'metadata'/('DOWNLOAD454_'+item['name']+'.json')
 if target.exists():
  assert target.stat().st_size==item['bytes'] and digest(target)==item['md5'];return dict(item,verified=True)
 for attempt in range(1,6):
  try:
   offset=part.stat().st_size if part.exists() else 0
   assert offset<=item['bytes']
   if offset<item['bytes']:
    headers={'Range':f'bytes={offset}-'} if offset else {}
    with session.get(item['url'],headers=headers,stream=True,timeout=(20,90)) as r:
     r.raise_for_status()
     if offset:assert r.status_code==206 and r.headers.get('Content-Range','').startswith(f'bytes {offset}-')
     with part.open('ab' if offset else 'wb') as f:
      last=0
      for chunk in r.iter_content(4194304):
       if chunk:f.write(chunk);offset+=len(chunk)
       assert offset<=item['bytes']
       if time.monotonic()-last>30:
        progress.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':offset,'total':item['bytes'],'attempt':attempt,'pid':os.getpid()}),encoding='utf-8');last=time.monotonic()
   assert part.stat().st_size==item['bytes']
   h=digest(part);assert h==item['md5'], 'MD5 mismatch; retain part for inspection'
   part.rename(target);return dict(item,verified=True,local_path=str(target))
  except Exception as e:
   print(item['name'],attempt,type(e).__name__,str(e),flush=True)
   if attempt==5:raise
   time.sleep(5)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:receipts=list(pool.map(download,m['items']))
(b/'metadata/PILOT_DOWNLOAD_COMPLETE_454.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'run':m['run'],'files':receipts,'all_md5_verified':True},indent=2),encoding='utf-8')
print('DOWNLOAD_COMPLETE',flush=True)
