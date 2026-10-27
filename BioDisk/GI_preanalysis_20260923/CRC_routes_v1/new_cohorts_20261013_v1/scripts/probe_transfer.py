from pathlib import Path
import requests,json,concurrent.futures,time
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1')
name='GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz'
urls=['http://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl/'+name,
'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200997&format=file&file='+name,
'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl/'+name+'?download=1']
def probe(url):
 start=time.time()
 try:
  s=requests.Session();s.trust_env=False
  with s.get(url,headers={'Range':'bytes=0-131071'},timeout=(10,15),stream=True) as r:
   r.raise_for_status();block=next(r.iter_content(16384))
   return dict(url=url,status=r.status_code,elapsed=time.time()-start,content_range=r.headers.get('Content-Range'),content_length=r.headers.get('Content-Length'),bytes=len(block),gzip=block[:2]==b'\x1f\x8b')
 except Exception as e:return dict(url=url,elapsed=time.time()-start,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:out=list(pool.map(probe,urls))
(R/'metadata/TRANSFER_ROUTE_PROBE.json').write_text(json.dumps(out,indent=2))
for x in out:print(json.dumps(x),flush=True)
