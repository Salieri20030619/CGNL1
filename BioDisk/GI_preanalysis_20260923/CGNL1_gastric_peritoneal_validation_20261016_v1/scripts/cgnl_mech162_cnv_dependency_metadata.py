from pathlib import Path
import requests,json,importlib.metadata,datetime
O=Path.cwd()/'mechanism_investigation_v2';u='https://pypi.org/pypi/infercnvpy/json';a=requests.get(u,timeout=(15,45));a.raise_for_status();q=a.json();out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':u,'package':q['info']['name'],'version':q['info']['version'],'requires_python':q['info']['requires_python'],'requires_dist':q['info']['requires_dist'],'files':[{'name':f['filename'],'url':f['url'],'sha256':f['digests']['sha256'],'size':f['size']} for f in q['urls']],'installed':{}}
for name in ['scanpy','anndata','numpy','pandas','scipy','joblib','matplotlib','tqdm','intervaltree','pyreadr','pybiomart','legacy-api-wrap']:
 try:out['installed'][name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:out['installed'][name]=None
(O/'metadata/CNV_DEPENDENCIES_162.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out));md=(O/'reference/gencode_v47/MD5SUMS').read_text(encoding='utf-8');print('GTF_CHECKSUMS',*[l for l in md.splitlines() if '.gtf.gz' in l and ('basic.annotation' in l or 'primary_assembly.annotation' in l)],sep='\n')
