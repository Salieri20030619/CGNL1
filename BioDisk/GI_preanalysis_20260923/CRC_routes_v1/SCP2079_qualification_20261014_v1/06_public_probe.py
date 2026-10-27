from pathlib import Path
import requests,json,concurrent.futures,hashlib
R=Path(__file__).resolve().parent
bucket=json.loads((R/'sources/explore.txt').read_text(encoding='utf-8'))['bucketId']
files=['expression_epi_counts/barcodes.tsv.gz','expression_epi_counts/matrix.mtx.gz']
def probe(name):
 url='https://storage.googleapis.com/'+bucket+'/'+name
 with requests.get(url,headers={'Range':'bytes=0-65535'},stream=True,timeout=30) as r:
  data=next(r.iter_content(65536),b'')
  d={'url':url,'status':r.status_code,'content_range':r.headers.get('Content-Range'),'content_type':r.headers.get('Content-Type'),'probe_bytes':len(data)}
  if r.status_code in (200,206):(R/'sources'/('public_probe_'+name.split('/')[-1])).write_bytes(data)
  print(json.dumps(d),flush=True);return d
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex: rows=list(ex.map(probe,files))
(R/'PUBLIC_OBJECT_PROBE.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')