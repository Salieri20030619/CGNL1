from pathlib import Path
import requests,json,hashlib,concurrent.futures,datetime,csv
b=Path(__file__).resolve().parents[1];urls=json.loads((b/'metadata/CPTAC_SOURCE_MANIFEST_349.json').read_text())['urls']
def fetch(u):
 name=u.rsplit('/',1)[1];p=b/'raw'/name;partial=p.with_suffix(p.suffix+'.partial');rec={'url':u,'name':name,'ok':False}
 try:
  with requests.get(u,stream=True,timeout=(15,60)) as r:
   r.raise_for_status();expected=int(r.headers.get('Content-Length',0));sha=hashlib.sha256();n=0
   with partial.open('wb') as f:
    for chunk in r.iter_content(1024*1024):f.write(chunk);sha.update(chunk);n+=len(chunk)
  assert not expected or n==expected
  partial.rename(p);rec.update(ok=True,bytes=n,sha256=sha.hexdigest())
  if any(k in name for k in ['proteomics','phospho']):
   hits=[]
   with p.open(encoding='utf-8') as f:
    header=f.readline().rstrip('\n');rowcount=0
    for line in f:
     rowcount+=1
     if 'CGNL1' in line.split('\t',1)[0]:hits.append(line.rstrip('\n'))
   # Qualification rows only; no contrasts or hypothesis tests executed
   (b/'metadata'/(name+'.CGNL1_rows_350.json')).write_text(json.dumps({'header':header,'total_rows':rowcount,'target_rows':hits},indent=2),encoding='utf-8');rec.update(total_rows=rowcount,CGNL1_rows=len(hits))
 except Exception as e:rec['error']=repr(e)
 print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(fetch,urls))
(b/'metadata/CPTAC_DOWNLOAD_COMPLETE_350.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
if not all(r['ok'] for r in out):raise SystemExit(1)
