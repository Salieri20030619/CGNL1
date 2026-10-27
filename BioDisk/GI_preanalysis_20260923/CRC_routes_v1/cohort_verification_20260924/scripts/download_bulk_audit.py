from pathlib import Path
import urllib.request,gzip,hashlib,json,time,traceback
D=Path(__file__).resolve().parent.parent/'bulk';D.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897','http':'http://127.0.0.1:7897'}))
urls={
'GSE159216_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159216/matrix/GSE159216_series_matrix.txt.gz',
'GSE73255_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73255/matrix/GSE73255_series_matrix.txt.gz',
'GPL17586_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17586/soft/GPL17586_family.soft.gz',
'GPL10558.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz'}
out=[]
for name,url in urls.items():
 p=D/name
 try:
  if not p.exists():
   for attempt in range(3):
    try:
     with op.open(url,timeout=90) as r,open(str(p)+'.part','wb') as f:
      size=r.headers.get('Content-Length');n=0
      while True:
       b=r.read(1024*1024)
       if not b:break
       f.write(b);n+=len(b)
     if size:assert n==int(size)
     Path(str(p)+'.part').replace(p);break
    except Exception:
     if attempt==2:raise
     time.sleep(2)
  with gzip.open(p,'rb') as f:
   while f.read(1024*1024):pass
  out.append({'file':name,'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'gzip_verified':True})
 except Exception:out.append({'file':name,'url':url,'error':traceback.format_exc()})
 (D/'MANIFEST.json').write_text(json.dumps(out,indent=2))
 print(name,'verified' if 'error' not in out[-1] else 'failed',flush=True)
