from pathlib import Path
import requests,gzip,json,concurrent.futures,threading,hashlib,datetime
O=Path.cwd()/'mechanism_investigation_v2';R=O/'raw/GSE237876';tasks=json.loads((O/'metadata/ROUTE_INPUT_PLAN_43.json').read_text());out=[];lock=threading.Lock()
def get(t):
 u=t['files'][0].replace('ftp:','http:');p=R/u.rsplit('/',1)[1];r=dict(t,url=u)
 try:
  if p.exists():b=p.read_bytes();gzip.decompress(b);r['reused_validated_existing']=True
  else:
   a=requests.get(u,timeout=(15,90));a.raise_for_status();b=a.content
   try:gzip.decompress(b).decode('utf-8')
   except Exception:(R/(p.name+'.HTTPfailed')).write_bytes(b);raise
   p.write_bytes(b)
  r.update(success=True,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:r.update(success=False,error=repr(e))
 with lock:
  out.append(r);(O/'metadata/ROUTE_HTTP_RECEIPTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(t['accession'],r['success'],r.get('error','CRCverified'),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rs=list(ex.map(get,tasks))
(O/'metadata/ROUTE_HTTP_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':all(r['success'] for r in rs),'passed':sum(r['success'] for r in rs),'files':len(rs)},indent=2),encoding='utf-8')
