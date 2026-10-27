from pathlib import Path
import requests,json,time,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';t=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text())['tasks'][0]
urls={'ENA_http':t['url'].replace('https:','http:'),'GENCODE_http':'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz'}
def probe(kv):
 k,u=kv;r={'key':k,'url':u};start=time.monotonic()
 try:
  with requests.get(u,headers={'Range':'bytes=0-1048575','Accept-Encoding':'identity'},stream=True,timeout=(15,45)) as a:
   b=a.raw.read(1048576);r.update(status=a.status_code,range=a.headers.get('Content-Range'),bytes=len(b),seconds=time.monotonic()-start,sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:r['error']=repr(e)
 print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(probe,urls.items()))
(O/'metadata/HTTP_CAPACITY_PROBE.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
