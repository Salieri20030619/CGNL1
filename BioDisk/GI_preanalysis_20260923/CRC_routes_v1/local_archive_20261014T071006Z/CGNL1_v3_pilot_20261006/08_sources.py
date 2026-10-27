from pathlib import Path
import requests,concurrent.futures,json,base64,hashlib,time
b=Path(r'C:\Users\31937\CGNL1_v3_pilot_20261006\direct_local')
tree=json.loads((b/'quiche_tree.json').read_text(encoding='utf-8'))
tasks=[(x['path'],x['url'],'blob') for x in tree['tree'] if x['path'].startswith('src/quiche/') and x['path'].endswith('.py')]
tasks += [('pyg_cpu.html','https://data.pyg.org/whl/torch-2.12.0+cpu.html','raw'),('pyg_cuda.html','https://data.pyg.org/whl/torch-2.12.0+cu130.html','raw'),('tissuumaps_tree.json','https://api.github.com/repos/TissUUmaps/TissUUmaps/git/trees/main?recursive=1','raw')]
def f(it):
 n,u,kind=it;s=requests.Session();s.trust_env=False;r={'name':n,'url':u,'network':'direct_only'}
 try:
  q=s.get(u,timeout=(8,25));q.raise_for_status();data=base64.b64decode(q.json()['content']) if kind=='blob' else q.content
  p=b/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data);r.update(status='downloaded',bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
 except Exception as e:r.update(status='failed',error=str(e))
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rs=list(ex.map(f,tasks))
(b/'source_fetch_log.json').write_text(json.dumps(rs,indent=2),encoding='utf-8');print(json.dumps(rs))
