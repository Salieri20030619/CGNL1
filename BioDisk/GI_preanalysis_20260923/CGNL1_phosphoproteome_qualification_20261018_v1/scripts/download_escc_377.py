from pathlib import Path
import requests,json,hashlib,zipfile,concurrent.futures
b=Path(__file__).resolve().parents[1];urls=json.loads((b/'metadata/ESCC_DOWNLOAD_MANIFEST_377.json').read_text())
def run(u):
 p=b/'raw'/u.rsplit('/',1)[1];sha=hashlib.sha256();n=0
 with requests.get(u,stream=True,timeout=(15,60)) as r:
  r.raise_for_status();expected=int(r.headers.get('Content-Length',0))
  with p.with_suffix('.partial').open('wb') as f:
   for chunk in r.iter_content(1048576):f.write(chunk);sha.update(chunk);n+=len(chunk)
 assert not expected or n==expected;p.with_suffix('.partial').rename(p)
 with zipfile.ZipFile(p) as z:assert z.testzip() is None
 return {'url':u,'path':str(p),'bytes':n,'sha256':sha.hexdigest(),'CRC_pass':True}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:out=list(ex.map(run,urls))
(b/'metadata/ESCC_DOWNLOAD_COMPLETE_377.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
