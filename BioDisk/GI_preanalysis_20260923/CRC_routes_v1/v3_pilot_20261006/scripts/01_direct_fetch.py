from pathlib import Path
import os,json,hashlib,datetime,time,concurrent.futures,subprocess,sys
ROOT=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006')
BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ['NO_PROXY']='*';os.environ['no_proxy']='*'
os.environ['PIP_CONFIG_FILE']='/dev/null'
import requests
tasks={
'zenodo_hd.json':'https://zenodo.org/api/records/19043057',
'iss_panel_2.tmap':'https://crclm2.serve.scilifelab.se/iss_panel_2.tmap',
'quiche_pyproject.toml':'https://raw.githubusercontent.com/jranek/quiche/main/pyproject.toml',
'quiche_tree.json':'https://api.github.com/repos/jranek/quiche/git/trees/main?recursive=1',
'nichecompass_pyproject.toml':'https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/main/pyproject.toml',
'nichecompass_tree.json':'https://api.github.com/repos/Lotfollahi-lab/nichecompass/git/trees/main?recursive=1',
'quiche_main.zip':'https://codeload.github.com/jranek/quiche/zip/refs/heads/main',
}
def get(it):
 name,url=it;s=requests.Session();s.trust_env=False;t=time.time();out=ROOT/'metadata'/name;part=out.with_suffix(out.suffix+'.part')
 r={'name':name,'url':url,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'network_mode':'direct_only','trust_env':False}
 try:
  resp=s.get(url,timeout=(12,35),stream=True);r.update(http_status=resp.status_code,final_url=resp.url);resp.raise_for_status()
  size=0
  with part.open('wb') as f:
   for ch in resp.iter_content(1024*1024):
    size+=len(ch)
    if size>100*1024*1024:raise ValueError('metadata/source cap exceeded')
    f.write(ch)
  declared=resp.headers.get('Content-Length')
  if declared and not resp.headers.get('Content-Encoding') and int(declared)!=size:raise ValueError('length mismatch')
  part.replace(out);r.update(status='downloaded',bytes=size,sha256=hashlib.sha256(out.read_bytes()).hexdigest())
 except Exception as e:r.update(status='failed',error=str(e))
 r['elapsed_s']=round(time.time()-t,3)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(get,tasks.items()))
(ROOT/'metadata/DIRECT_FETCH_LOG.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
(ROOT/'NETWORK_MODE.json').write_text(json.dumps({'mode':'direct_only','requests_trust_env':False,'proxy_variables_removed':True,'NO_PROXY':'*','fallback_to_proxy':False},indent=2))
for r in results:print(json.dumps(r),flush=True)
m=json.loads((ROOT/'PLAN_FREEZE.json').read_text())
assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in m['files'].items())
(ROOT/'metadata/PLAN_FREEZE_VERIFIED.json').write_text(json.dumps({'passed':True,'checked_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}))
print('DONE_DIRECT_FETCH',flush=True)
