import pathlib, urllib.request, gzip, json, hashlib, time, traceback
D=pathlib.Path(__file__).parent/'reference_supplement_20260924';D.mkdir(exist_ok=True)
S=D/'STATE.json'
def state(**kw):S.write_text(json.dumps(kw,indent=2,ensure_ascii=False),encoding='utf-8')
op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}))
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/'
files={'GSE115469_CellClusterType.txt.gz':base+'suppl/GSE115469_CellClusterType.txt.gz','GSE115469_family.soft.gz':base+'soft/GSE115469_family.soft.gz','GSE115469_Data.csv.gz':base+'suppl/GSE115469_Data.csv.gz'}
manifest=[]
try:
 for name,url in files.items():
  p=D/name;state(status='downloading',file=name,completed=manifest)
  for attempt in range(4):
   try:
    if not p.exists():
     with op.open(url,timeout=60) as r,open(str(p)+'.part','wb') as w:
      expected=r.headers.get('Content-Length');size=0
      while True:
       b=r.read(1024*1024)
       if not b:break
       w.write(b);size+=len(b)
     if expected and size!=int(expected):raise ValueError('Content length mismatch')
     pathlib.Path(str(p)+'.part').replace(p)
    with gzip.open(p,'rb') as f:
     while f.read(1024*1024):pass
    break
   except Exception:
    if attempt==3:raise
    time.sleep(3*(attempt+1))
  h=hashlib.sha256()
  with open(p,'rb') as f:
   while True:
    b=f.read(1024*1024)
    if not b:break
    h.update(b)
  manifest.append({'file':name,'url':url,'bytes':p.stat().st_size,'sha256':h.hexdigest(),'gzip_verified':True})
  with gzip.open(p,'rt',encoding='utf-8') as f:preview=''.join(f.readline() for _ in range(4))
  (D/(name+'.preview.txt')).write_text(preview[:50000],encoding='utf-8')
 (D/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 state(status='download_verified_pending_reference_review',completed=manifest,next_steps=['Verify author hepatocyte annotation and donor mapping','Verify expression scale; generate equal-donor reference profiles','Keep >=3 donor and >=50 cell profile gate','Resume composition and external spatial stages in new version; preserve previous results'],scope='Do not merge into analysis until donor mapping and expression scale are reviewed')
 (D/'EXIT.json').write_text(json.dumps({'exit_code':0,'completion':'download_and_gzip_integrity_verified'}))
except Exception:
 state(status='failed',traceback=traceback.format_exc(),completed=manifest)
 (D/'EXIT.json').write_text(json.dumps({'exit_code':1}))
 raise
