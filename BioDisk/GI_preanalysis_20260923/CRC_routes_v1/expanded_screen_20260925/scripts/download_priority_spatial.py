from pathlib import Path
import urllib.request,concurrent.futures,gzip,json,hashlib
A=Path(__file__).resolve().parent.parent;M=A/'metadata';D=A/'spatial_data';D.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}));jobs=[]
for g in ['GSE206552','GSE254364','GSE217414']:
 p=M/(g+'.soft.gz')
 if not p.exists():
  with op.open('https://ftp.ncbi.nlm.nih.gov/geo/series/'+g[:-3]+'nnn/'+g+'/soft/'+g+'_family.soft.gz',timeout=45) as f:p.write_bytes(f.read())
 lines=gzip.decompress(p.read_bytes()).decode().splitlines();s=None;meta=[]
 for l in lines:
  if l.startswith('^SAMPLE'):s={'gsm':l.split(' = ')[1],'urls':[]};meta.append(s)
  if s is not None:
   if l.startswith('!Sample_title'):s['title']=l.split(' = ')[1]
   if l.startswith('!Sample_characteristics'):s.setdefault('characteristics',[]).append(l.split(' = ')[1])
   if l.startswith('!Sample_supplementary_file'):
    u=l.split(' = ')[1].replace('ftp://','https://');s['urls'].append(u)
    if u!='NONE' and 'raw_feature' not in u:jobs.append((g,u))
 (M/(g+'_samples.json')).write_text(json.dumps(meta,indent=2))
 print(g,len(meta),'samples',sum(len(x['urls']) for x in meta),'files',flush=True)
def download(j):
 g,u=j;p=D/g/u.rsplit('/',1)[1];p.parent.mkdir(exist_ok=True)
 try:
  if not p.exists():
   with op.open(u,timeout=60) as r,open(str(p)+'.part','wb') as f:
    size=r.headers.get('Content-Length');n=0
    if size:assert int(size)<300000000
    while True:
     b=r.read(1048576)
     if not b:break
     f.write(b);n+=len(b);assert n<300000000
   if size:assert int(size)==n
   Path(str(p)+'.part').replace(p)
  if p.suffix=='.gz':
   with gzip.open(p,'rb') as f:
    while f.read(1048576):pass
  rec={'dataset':g,'url':u,'file':str(p.relative_to(A)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'status':'verified'};print(p.name,'OK',flush=True);return rec
 except Exception as e:print(u,str(e),flush=True);return {'dataset':g,'url':u,'status':'error','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(download,jobs))
(A/'SPATIAL_DOWNLOAD_MANIFEST.json').write_text(json.dumps(res,indent=2));(A/'SPATIAL_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'exit_code':0 if all(x['status']=='verified' for x in res) else 1,'files':len(res)}));print('DOWNLOAD_COMPLETE',flush=True)
