import requests,time,json,concurrent.futures
from pathlib import Path
L=Path(__file__).resolve().parent
url='http://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz'
def test(proxy):
 start=time.time();s=requests.Session();s.trust_env=False
 if proxy:s.proxies={'http':'http://127.0.0.1:7897'}
 try:
  with s.get(url,headers={'Range':'bytes=0-131071'},stream=True,timeout=(8,15)) as r:
   r.raise_for_status();b=r.raw.read(16384)
   return dict(proxy=proxy,status=r.status_code,bytes=len(b),elapsed=time.time()-start,gzip=b[:2]==b'\x1f\x8b',content_range=r.headers.get('Content-Range'))
 except Exception as e:return dict(proxy=proxy,elapsed=time.time()-start,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:res=list(pool.map(test,[False,True]))
(L/'metadata/LOCAL_HTTP_PROBE.json').write_text(json.dumps(res,indent=2))
for x in res:print(json.dumps(x),flush=True)
