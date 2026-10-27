from pathlib import Path
import gzip,re,json,urllib.request,hashlib,time,concurrent.futures,threading,traceback
D=Path(__file__).resolve().parent.parent;R=D/'GSE245552';R.mkdir(exist_ok=True)
s=gzip.open(D/'sc_metadata/GSE245552.soft.gz','rt',encoding='utf-8').read();samples=[]
for block in s.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];get=lambda prefix:[x.split(' = ',1)[1] for x in block.splitlines() if x.startswith(prefix)]
 title=get('!Sample_title')[0];tissue=get('!Sample_source_name_ch1')[0]
 if tissue!='liver metastasis':continue
 urls=get('!Sample_supplementary_file_');urls=[u.replace('ftp://','https://') for u in urls if u!='NONE'];assert len(urls)==3
 donor=re.match(r'(Pt\d+)_',title).group(1);samples.append({'dataset':'GSE245552','sample':gsm,'donor':donor,'title':title,'tissue':tissue,'urls':urls})
(D/'GSE245552_LM_SAMPLES.json').write_text(json.dumps(samples,indent=2))
print('Liver metastasis libraries:',len(samples),'donor labels:',len(set(x['donor'] for x in samples)),flush=True)
jobs=[u for x in samples for u in x['urls']]
def fetch(url):
 p=R/url.rsplit('/',1)[1];op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
 for attempt in range(3):
  try:
   if not p.exists():
    with op.open(url,timeout=90) as r,open(str(p)+'.part','wb') as f:
     size=r.headers.get('Content-Length');n=0
     while True:
      b=r.read(1024*1024)
      if not b:break
      f.write(b);n+=len(b)
    if size:assert n==int(size)
    Path(str(p)+'.part').replace(p)
   with gzip.open(p,'rb') as f:
    while f.read(1024*1024):pass
   return {'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'url':url,'gzip_verified':True}
  except Exception:
   if attempt==2:return {'file':p.name,'url':url,'error':traceback.format_exc()}
   time.sleep(2)
out=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for item in ex.map(fetch,jobs):
  out.append(item);(D/'GSE245552_DOWNLOAD_STATUS.json').write_text(json.dumps({'completed':len(out),'total':len(jobs),'files':out},indent=2));print(item['file'],'ok' if 'error' not in item else 'failed',flush=True)
assert all('error' not in x for x in out)
(D/'GSE245552_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'exit_code':0,'libraries':len(samples),'donor_labels':len(set(x['donor'] for x in samples)),'files':len(out),'bytes':sum(x['bytes'] for x in out)},indent=2))
