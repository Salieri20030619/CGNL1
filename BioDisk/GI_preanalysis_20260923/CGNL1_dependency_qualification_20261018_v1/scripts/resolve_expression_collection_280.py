from pathlib import Path
import json,hashlib,requests,threading,concurrent.futures,time,datetime,os
p=Path.cwd();original=json.loads((p/'metadata/TARGET_RNA_RECEIPTS_269.json').read_text());assert len(original)==162
recovery=json.loads((p/'metadata/RNA_RECOVERY_RECEIPT_279.json').read_text());out=p/'raw/target_rnaseq_resolved_280';out.mkdir(exist_ok=False);tls=threading.local()
params={'page[size]':0,'filter':json.dumps([{'name':'gene_id','op':'in_','val':['SIDG04503','SIDG04502','SIDG04102']}])}
def validate(raw,mid):
 d=json.loads(raw);assert isinstance(d['data'],list) and not d.get('links',{}).get('next')
 for row in d['data']:assert row['attributes']['model_id']==mid and row['attributes']['ensembl_gene_id'] in ['ENSG00000128849','ENSG00000143375','ENSG00000039068']
 return len(d['data'])
def resolve(r):
 mid=r['model_id'];dest=out/(mid+'.json');error=None
 if r['success']:
  source=p/'raw/target_rnaseq_269'/(mid+'.json');raw=source.read_bytes();assert hashlib.sha256(raw).hexdigest()==r['sha256'];n=validate(raw,mid);os.link(source,dest);return dict(model_id=mid,success=True,origin='original_verified',records=n,sha256=r['sha256'],source_path=str(source))
 if recovery['success'] and recovery['model_id']==mid:
  source=Path(recovery['file']);raw=source.read_bytes();n=validate(raw,mid);os.link(source,dest);return dict(model_id=mid,success=True,origin='separate_recovery279',records=n,sha256=hashlib.sha256(raw).hexdigest(),source_path=str(source))
 if not hasattr(tls,'session'):tls.session=requests.Session();tls.session.trust_env=False
 for attempt in range(1,3):
  try:
   response=tls.session.get('https://api.cellmodelpassports.sanger.ac.uk/models/'+mid+'/datasets/rnaseq',params=params,timeout=(15,60));response.raise_for_status();raw=response.content;n=validate(raw,mid);dest.write_bytes(raw);return dict(model_id=mid,success=True,origin='bounded_recovery280',records=n,sha256=hashlib.sha256(raw).hexdigest(),attempt=attempt)
  except Exception as exc:
   error=repr(exc)
   if isinstance(exc,requests.HTTPError) and exc.response.status_code in [401,403]:break
   if attempt<2:time.sleep(15)
 return dict(model_id=mid,success=False,error=error)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:receipts=list(ex.map(resolve,original))
(p/'metadata/RESOLVED_RNA_RECEIPTS_280.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8');x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':len(receipts),'all_requests_pass':all(r['success'] for r in receipts),'original_receipts_preserved':True,'recovered_models':[r['model_id'] for r in receipts if r.get('origin')!='original_verified']};(p/'metadata/RESOLVED_RNA_COMPLETE_280.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True);assert x['all_requests_pass']
