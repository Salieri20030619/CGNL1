from pathlib import Path
import requests,json,hashlib,time
O=Path.cwd()/'mechanism_investigation_v2';p=O/'raw/GSE308231_RAW.tar.partial';u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308231/suppl/GSE308231_RAW.tar';r={'url':u,'range_requested':'bytes=0-65535','local_prefix_bytes':p.stat().st_size if p.exists() else None};start=time.monotonic()
try:
 with requests.get(u,headers={'Range':'bytes=0-65535','Accept-Encoding':'identity'},stream=True,timeout=(15,40)) as a:
  b=a.raw.read(65536);r.update(status=a.status_code,content_range=a.headers.get('Content-Range'),response_bytes=len(b),elapsed=time.monotonic()-start,etag=a.headers.get('ETag'),last_modified=a.headers.get('Last-Modified'))
 with p.open('rb') as f:r['matches_existing_prefix']=b==f.read(len(b))
 r['safe_range_capability']=r['status']==206 and r['matches_existing_prefix'] and len(b)==65536
except Exception as e:r['error']=repr(e)
(O/'metadata/GSE308231_RANGE_PROBE.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
