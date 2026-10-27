from pathlib import Path
import requests,json,hashlib,concurrent.futures,time
b=Path(r'C:\Users\31937\CGNL1_v3_pilot_20261006\direct_local');w=b/'wheels';w.mkdir(exist_ok=True)
pkgs={'nichecompass':'0.3.3','liana':'1.10.0','torch-geometric':None,'mlflow-skinny':'2.22.2','pyreadr':None,'decoupler':'1.9.2','plotnine':'0.15.0','mizani':'0.14.2','xxhash':None,'session-info2':None,'quiche':None}
def run(it):
 n,v=it;s=requests.Session();s.trust_env=False;r={'package':n,'mode':'direct_only','trust_env':False}
 try:
  u='https://pypi.org/pypi/'+n+('/'+v if v else '')+'/json';j=s.get(u,timeout=(10,25));j.raise_for_status();d=j.json();(b/(n+'_metadata.json')).write_text(json.dumps(d),encoding='utf-8')
  r.update(version=d['info']['version'],project_urls=d['info'].get('project_urls'),requires_dist=d['info'].get('requires_dist'))
  if n=='quiche':return r
  options=[x for x in d['urls'] if x['filename'].endswith('.whl') and (x['filename'].endswith('py3-none-any.whl') or x['filename'].endswith('py2.py3-none-any.whl') or ('cp311' in x['filename'] and 'manylinux' in x['filename'] and 'x86_64' in x['filename']))]
  if not options:raise ValueError('no compatible Linux cp311 wheel')
  f=sorted(options,key=lambda x:x['size'])[0];q=s.get(f['url'],timeout=(10,40));q.raise_for_status();data=q.content
  assert hashlib.sha256(data).hexdigest()==f['digests']['sha256']
  (w/f['filename']).write_bytes(data);r.update(file=f['filename'],bytes=len(data),sha256=f['digests']['sha256'],url=f['url'],status='verified')
 except Exception as e:r.update(status='failed',error=str(e))
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:r=list(ex.map(run,pkgs.items()))
(b/'wheel_download_log.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
for x in r:print(json.dumps({k:v for k,v in x.items() if k not in ['requires_dist']}))
