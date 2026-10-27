from pathlib import Path
import urllib.request,json,hashlib
A=Path(__file__).resolve().parent.parent;D=A/'multiome';D.mkdir(exist_ok=True);plan=json.loads((A/'MULTIOME_DOWNLOAD_PLAN.json').read_text());op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}));out=[]
for x in plan:
 p=D/x['filename'];info=x['content_details']
 if not p.exists():
  with op.open(info['download_url'],timeout=90) as r,open(str(p)+'.part','wb') as f:
   while True:
    b=r.read(1048576)
    if not b:break
    f.write(b)
  Path(str(p)+'.part').replace(p)
 assert p.stat().st_size==info['size'];h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==info['sha256_hash'];out.append({'file':p.name,'bytes':p.stat().st_size,'sha256':h,'publisher_hash_verified':True});print(p.name,'verified',flush=True)
(A/'MULTIOME_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'exit_code':0,'files':out},indent=2))
