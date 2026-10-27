from pathlib import Path
import urllib.request,json,gzip,hashlib
D=Path(__file__).resolve().parent.parent/'GSE178318';D.mkdir(exist_ok=True);op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}));out=[]
for name in ['GSE178318_barcodes.tsv.gz','GSE178318_genes.tsv.gz','GSE178318_matrix.mtx.gz']:
 p=D/name;url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/'+name
 if not p.exists():
  with op.open(url,timeout=90) as r,open(str(p)+'.part','wb') as f:
   size=r.headers.get('Content-Length');n=0
   if size:assert int(size)<1500*1024*1024
   while True:
    b=r.read(1024*1024)
    if not b:break
    f.write(b);n+=len(b);assert n<1500*1024*1024
  if size:assert n==int(size)
  Path(str(p)+'.part').replace(p)
 with gzip.open(p,'rb') as f:
  while f.read(1024*1024):pass
 out.append({'file':name,'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'gzip_verified':True});(D/'MANIFEST.json').write_text(json.dumps(out,indent=2));print(name,'verified',flush=True)
with gzip.open(D/'GSE178318_barcodes.tsv.gz','rt') as f:print('barcode preview:',''.join(f.readline() for _ in range(4)),flush=True)
with gzip.open(D/'GSE178318_genes.tsv.gz','rt') as f:print('feature preview:',''.join(f.readline() for _ in range(4)),flush=True)
(D/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'files':3}))
