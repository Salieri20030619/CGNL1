from pathlib import Path
import json,urllib.request,gzip,hashlib,datetime,concurrent.futures,threading
O=Path.cwd()/'mechanism_investigation_v2';R=O/'raw/GSE237876';R.mkdir(exist_ok=True);plan=json.loads((O/'metadata/ROUTE_INPUT_PLAN_43.json').read_text());out=[];lock=threading.Lock()
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'route-specific tissue association, not cell-intrinsic mechanism','target':'CGNL1','units':'14patients,5withPM; collapse replicates and regions within patient/site before paired comparisons','primary':'patient-mean PM minus patient-mean primary expression; n5 descriptive exploratory','secondary':'other metastatic sites minus primary; report patient overlap and do not treat routes as independent when same donor','source_proxies':'fixed prior6lineage panels descriptive with low n; no overfit adjustment','qualification':'inspect actual scale and target identifiers before expression analysis; no missing-as-zero'}
(O/'ROUTE_ANALYSIS_INTENT.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
def get(t):
 u=t['files'][0];n=u.rsplit('/',1)[1];p=R/n;r=dict(t,url=u)
 try:
  assert not p.exists(),'existingfilemustinspect'
  with urllib.request.urlopen(u,timeout=60) as a:b=a.read(20*1024*1024+1)
  assert len(b)<=20*1024*1024
  try:raw=gzip.decompress(b);raw.decode('utf-8')
  except Exception:
   (R/(n+'.failed')).write_bytes(b);raise
  p.write_bytes(b);r.update(success=True,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),uncompressed_bytes=len(raw),header=raw.splitlines()[0].decode('utf-8'))
 except Exception as e:r.update(success=False,error=repr(e))
 with lock:
  out.append(r);(O/'metadata/ROUTE_DOWNLOAD_RECEIPTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(t['accession'],r['success'],r.get('header',r.get('error')),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rs=list(ex.map(get,plan))
(O/'metadata/ROUTE_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'files':len(rs),'success':sum(r['success'] for r in rs),'all_pass':all(r['success'] for r in rs)},indent=2),encoding='utf-8')
