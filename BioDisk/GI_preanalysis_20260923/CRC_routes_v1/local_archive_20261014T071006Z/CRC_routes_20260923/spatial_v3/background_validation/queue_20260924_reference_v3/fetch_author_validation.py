from pathlib import Path
import urllib.request,json,hashlib,traceback
D=Path(__file__).parent/'AUTHOR_VALIDATION';D.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897','http':'http://127.0.0.1:7897'}))
items=[('integrated_metCRC_downstream_analysis.R','https://ndownloader.figshare.com/files/43628790','357c9476013f3a5208d571e42109191d'),('README.md','https://ndownloader.figshare.com/files/43628793','c03ddaa19d0f85d4fbcab1bad73ae1f3'),('merged_NDE_22261_fisher.rds','https://ndownloader.figshare.com/files/43628796','149423456cc48658d2cdb51a8908d2b8'),('author_spatial_deconvolution.pdf','https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13059-023-03159-6/MediaObjects/13059_2023_3159_MOESM3_ESM.pdf',None)]
rows=[]
for name,url,md5 in items:
 p=D/name
 try:
  if not p.exists():
   with op.open(url,timeout=90) as r,open(str(p)+'.part','wb') as f:
    while True:
     b=r.read(1024*1024)
     if not b:break
     f.write(b)
   Path(str(p)+'.part').replace(p)
  if md5:assert hashlib.md5(p.read_bytes()).hexdigest()==md5
  rows.append({'file':name,'bytes':p.stat().st_size,'url':url,'md5':hashlib.md5(p.read_bytes()).hexdigest(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'status':'verified' if md5 else 'downloaded'})
 except Exception:rows.append({'file':name,'status':'failed','traceback':traceback.format_exc()})
 (D/'MANIFEST.json').write_text(json.dumps(rows,indent=2))
 print(name,rows[-1]['status'],flush=True)
