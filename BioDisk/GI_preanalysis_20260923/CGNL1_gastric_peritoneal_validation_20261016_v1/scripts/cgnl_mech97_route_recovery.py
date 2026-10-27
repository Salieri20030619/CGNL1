from pathlib import Path
import requests,json,gzip,hashlib,datetime,concurrent.futures,threading,os,signal
O=Path.cwd()/'mechanism_investigation_v2';assert not (O/'results/ROUTE_ANALYSIS_COMPLETE.json').exists();old=json.loads((O/'metadata/ROUTE_HTTP_RECEIPTS.json').read_text());failed=[r for r in old if not r['success']];assert len(failed)==20;stopped=[]
for p in Path('/proc').iterdir():
 if p.name.isdigit():
  try:
   if str(O.parent/'scripts/cgnl_mech71_route_analyze_when_ready.py').encode() in (p/'cmdline').read_bytes().split(b'\0'):os.kill(int(p.name),signal.SIGTERM);stopped.append(int(p.name))
  except (FileNotFoundError,PermissionError):pass
(O/'metadata/ROUTE_RECOVERY_PLAN.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'initial20HTTP503failuresfollowedbycontinuoussuccessfuldownloads;oneboundedretryofonlyfailedfilesbeforeexpressionoutcomes','stopped_pending_analysis':stopped,'failed_samples':[r['accession'] for r in failed],'prior_receipts':old},indent=2),encoding='utf-8');out=[];lock=threading.Lock()
def get(t):
 r=dict(t);u=t['files'][0].replace('ftp:','http:');p=O/'raw/GSE237876'/u.rsplit('/',1)[1]
 try:
  assert not p.exists(),'existingfilemustinspect';a=requests.get(u,timeout=(15,90));a.raise_for_status();b=a.content;gzip.decompress(b).decode('utf-8');p.write_bytes(b);r.update(success=True,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest());r.pop('error',None)
 except Exception as e:r.update(success=False,error=repr(e))
 with lock:out.append(r);(O/'metadata/ROUTE_RECOVERY_RECEIPTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(t['accession'],r['success'],r.get('error','CRCverified'),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:new=list(ex.map(get,failed))
combined={r['accession']:r for r in old};combined.update({r['accession']:r for r in new});combined=list(combined.values());(O/'metadata/ROUTE_HTTP_FINAL_RECEIPTS.json').write_text(json.dumps(combined,indent=2),encoding='utf-8');(O/'metadata/ROUTE_HTTP_RECOVERY_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':66,'passed':sum(r['success'] for r in combined),'all_pass':all(r['success'] for r in combined),'one_recovery_round':True},indent=2),encoding='utf-8')
