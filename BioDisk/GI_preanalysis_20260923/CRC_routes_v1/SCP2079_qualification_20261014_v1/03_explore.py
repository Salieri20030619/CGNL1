from pathlib import Path
import requests,json,concurrent.futures
ROOT=Path(__file__).resolve().parent
base='https://singlecell.broadinstitute.org/single_cell/api/v1/'
urls={'explore':base+'studies/SCP2079/explore','annotations':base+'studies/SCP2079/annotations','clusters':base+'studies/SCP2079/clusters','download_probe':base+'site/studies/SCP2079/download?filename=epi_metadata.txt'}
def f(k,u):
 try:
  r=requests.get(u,timeout=30); (ROOT/'sources'/f'{k}.txt').write_bytes(r.content)
  print(k,r.status_code,len(r.content),flush=True)
  if k!='explore':print(r.text[:6000],flush=True)
  else:
   d=r.json();print('Explore keys:',list(d),flush=True)
   print(json.dumps({k:v for k,v in d.items() if not any(w in k.lower() for w in ['token','auth'])})[:10000],flush=True)
  return dict(name=k,url=u,status=r.status_code,bytes=len(r.content))
 except Exception as e:return dict(name=k,error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: rows=list(ex.map(lambda x:f(*x),urls.items()))
(ROOT/'EXPLORE_FETCH.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')