from pathlib import Path
import urllib.request,gzip,re,json,concurrent.futures,time
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');(O/'raw').mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
jobs=[]
for g in ['GSE217414','GSE335552']:
 s=gzip.decompress((O/'metadata'/(g+'.soft.gz')).read_bytes()).decode()
 for url in re.findall(r'!Sample_supplementary_file(?:_\d+)? = (ftp[^\n]+\.h5)',s):
  url=url.replace('ftp://','https://');jobs.append((url,O/'raw'/url.rsplit('/',1)[-1]))
s=(O/'metadata/BJC2026.html').read_text(encoding='utf-8');u=re.search(r'https://media[^\s<>"\']+\.xlsx',s).group();jobs.append((u,O/'metadata/BJC2026_gene_panels.xlsx'))
def fetch(job):
 u,p=job
 try:
  if not p.exists():
   start=time.monotonic()
   with op.open(u,timeout=30) as r,Path(str(p)+'.part').open('wb') as f:
    while b:=r.read(1024*512):
     f.write(b)
     if time.monotonic()-start>150:raise TimeoutError('150s overall limit')
   Path(str(p)+'.part').replace(p)
  return {'file':p.name,'url':u,'bytes':p.stat().st_size,'ok':True}
 except Exception as e:return {'file':p.name,'url':u,'ok':False,'error':str(e)}
res=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(fetch,jobs));(O/'metadata/MATRIX_FETCH.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
