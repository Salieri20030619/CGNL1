from pathlib import Path
import requests,concurrent.futures,json,time
b=Path(r'C:\Users\31937\CGNL1_v3_pilot_20261006');out=b/'direct_local';out.mkdir(exist_ok=True)
base='https://crclm2.serve.scilifelab.se/'
urls={
'iss_download':'download/iss_panel_2/prepared_cells_panel_2.h5ad',
'iss_dl':'iss_panel_2/prepared_cells_panel_2.h5ad?dl=1',
'iss_js':'static/js/utils/h5Utils.js',
'tissuumaps_views':'https://raw.githubusercontent.com/TissUUmaps/TissUUmaps/master/tissuumaps/views.py',
'pypi_nichecompass':'https://pypi.org/pypi/nichecompass/0.3.3/json',
'pypi_liana':'https://pypi.org/pypi/liana/json',
'pypi_torch_geo':'https://pypi.org/pypi/torch-geometric/json',
}
def f(it):
 n,p=it;u=p if p.startswith('https') else base+p;s=requests.Session();s.trust_env=False;t=time.time()
 r={'name':n,'url':u,'trust_env':False,'mode':'direct_only'}
 try:
  q=s.get(u,timeout=(8,20),stream=True);r.update(status=q.status_code,content_type=q.headers.get('Content-Type'),length=q.headers.get('Content-Length'));q.raise_for_status()
  if n in ['iss_download','iss_dl']:
   chunk=next(q.iter_content(64));r['magic_hex']=chunk[:16].hex();q.close()
  else:
   data=q.content;(out/(n+'.txt')).write_bytes(data);r['bytes']=len(data)
 except Exception as e:r['error']=str(e)
 r['elapsed']=time.time()-t;return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rs=list(ex.map(f,urls.items()))
(out/'probe_log.json').write_text(json.dumps(rs,indent=2));print(json.dumps(rs))
