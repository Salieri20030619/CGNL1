import requests,json,concurrent.futures
from pathlib import Path
L=Path(__file__).resolve().parent
urls={'khaliq_zenodo':'https://zenodo.org/api/records/6466249','lee_zenodo':'https://zenodo.org/api/records/3967538','khaliq_tree':'https://api.github.com/repos/masoodlab/CRC-Single-Cell/git/trees/0.1.1?recursive=1'}
def f(item):
 k,u=item
 try:
  r=requests.get(u,timeout=(10,25));r.raise_for_status();(L/'metadata'/(k+'.json')).write_text(r.text);j=r.json()
  if 'files' in j:return {k:[{'key':x['key'],'size':x['size'],'url':x['links']['self']} for x in j['files']]}
  return {k:[{'path':x['path'],'size':x.get('size',0)} for x in j.get('tree',[]) if x.get('size',0)>100000 or any(t in x['path'].lower() for t in ['matrix','seurat','metadata','annotation'])]}
 except Exception as e:return {k:repr(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for x in pool.map(f,urls.items()):print(json.dumps(x),flush=True)