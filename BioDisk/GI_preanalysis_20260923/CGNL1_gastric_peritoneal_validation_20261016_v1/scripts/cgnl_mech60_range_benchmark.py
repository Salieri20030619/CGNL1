from pathlib import Path
import requests,json,time,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';t=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text())['tasks'][0]
def test(n):
 r={'bytes_requested':n};start=time.monotonic()
 try:
  with requests.get(t['url'],headers={'Range':f'bytes=0-{n-1}','Accept-Encoding':'identity'},timeout=(15,45)) as a:
   b=a.content;r.update(status=a.status_code,content_range=a.headers.get('Content-Range'),bytes_received=len(b),elapsed=time.monotonic()-start);assert a.status_code==206 and len(b)==n
  p=O/'raw/DRA014736'/(t['name']+'.partial')
  with p.open('rb') as f:prefix=f.read(n)
  r['matches_existing_https_prefix']=prefix==b;r['sha256']=hashlib.sha256(b).hexdigest()
 except Exception as e:r['error']=repr(e)
 print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(test,[1024*1024,4*1024*1024]))
(O/'metadata/ENA_RANGE_BENCHMARK.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
