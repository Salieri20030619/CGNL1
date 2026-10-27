from pathlib import Path
import requests,json,hashlib,concurrent.futures
b=Path(__file__).resolve().parents[1];d=b/'raw/GC2026_410';d.mkdir(exist_ok=True);urls=json.loads((b/'metadata/GC2026_DOWNLOAD_MANIFEST_410.json').read_text())
def run(u):
 p=d/u.rsplit('/',1)[1];n=0;sha=hashlib.sha256()
 with requests.get(u,stream=True,timeout=(15,60)) as r:
  r.raise_for_status();expected=int(r.headers.get('Content-Length',0))
  with p.with_suffix(p.suffix+'.partial').open('wb') as f:
   for chunk in r.iter_content(1048576):f.write(chunk);sha.update(chunk);n+=len(chunk)
 assert not expected or n==expected;p.with_suffix(p.suffix+'.partial').rename(p);out={'url':u,'path':str(p),'bytes':n,'sha256':sha.hexdigest()};print(json.dumps(out),flush=True);return out
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(run,urls))
(b/'metadata/GC2026_DOWNLOAD_COMPLETE_410.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
