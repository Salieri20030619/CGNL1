from pathlib import Path
import requests,concurrent.futures,json,re,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1')
j=json.loads((R/'metadata/atlas_tree.json').read_text());sha=j['sha'];D=R/'metadata/atlas_loaders';D.mkdir(exist_ok=True)
paths=[x['path'] for x in j['tree'] if x['path'].startswith('analyses/01_dataloader/') and x['path'].endswith('.py')]
def fetch(path):
 url='https://raw.githubusercontent.com/icbi-lab/crc-atlas/'+sha+'/'+path;p=D/Path(path).name
 try:
  if not p.exists():
   s=requests.Session();s.trust_env=False
   r=s.get(url,timeout=(12,20));r.raise_for_status();p.write_bytes(r.content)
  txt=p.read_text();lines=txt.splitlines()
  hits=[dict(line=i+1,text=l[:240]) for i,l in enumerate(lines) if re.search(r'liver|metasta|mCRC|CRCLM',l,re.I)]
  return dict(path=path,url=url,status='read',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),accessions=sorted(set(re.findall(r'GSE\d+|E-MTAB-\d+|PRJNA\d+',txt))),liver_or_metastasis_mentions=hits)
 except Exception as e:return dict(path=path,url=url,status='failed',error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(fetch,paths))
(R/'metadata/ATLAS_SOURCE_AUDIT.json').write_text(json.dumps(dict(commit=sha,loaders=len(paths),records=res),indent=2))
for r in res:
 if r.get('liver_or_metastasis_mentions'):print(json.dumps(r),flush=True)
print('ATLAS_AUDIT_COMPLETE',len(res),flush=True)
