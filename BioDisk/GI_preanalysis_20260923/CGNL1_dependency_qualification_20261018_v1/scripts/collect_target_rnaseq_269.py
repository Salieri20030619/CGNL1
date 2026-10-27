from pathlib import Path
import json,csv,requests,concurrent.futures,time,datetime,hashlib
p=Path.cwd();out=p/'raw/target_rnaseq_269';out.mkdir(exist_ok=True)
with (p/'results/organoid_264/target_model_values.tsv').open(encoding='utf-8') as f:models=[r['model_id'] for r in csv.DictReader(f,delimiter='\t')]
params={'page[size]':0,'filter':json.dumps([{'name':'gene','op':'has','val':{'name':'symbol','op':'in_','val':['CGNL1','CGN','CDH1']}}])}
def fetch(mid):
 u='https://api.cellmodelpassports.sanger.ac.uk/models/'+mid+'/datasets/rnaseq';s=requests.Session();s.trust_env=False
 for attempt in range(1,4):
  try:
   if mid=='SIDM01281':raw=(p/'metadata/RNASEQ_TARGET_API_PROBE_267.json').read_bytes();d=json.loads(raw)
   else:
    r=s.get(u,params=params,timeout=(15,45));r.raise_for_status();raw=r.content;d=r.json()
   assert isinstance(d['data'],list)
   for row in d['data']:assert row['attributes']['model_id']==mid and row['attributes']['ensembl_gene_id'] in ['ENSG00000128849','ENSG00000143375','ENSG00000039068']
   assert not d.get('links',{}).get('next')
   f=out/(mid+'.json');f.write_bytes(raw);return {'model_id':mid,'success':True,'records':len(d['data']),'url':u,'sha256':hashlib.sha256(raw).hexdigest()}
  except Exception as e:
   if isinstance(e,requests.HTTPError) and e.response.status_code in [401,403]:break
   if attempt<3:time.sleep(10*attempt)
 return {'model_id':mid,'success':False,'error':repr(e) if 'e' in locals() else 'bounded request attempts exhausted'}
receipts=[]
def save():
 f=p/'metadata/TARGET_RNA_RECEIPTS_269.json';tmp=f.with_suffix('.writing');tmp.write_text(json.dumps(receipts,indent=2),encoding='utf-8');tmp.replace(f)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for future in concurrent.futures.as_completed([ex.submit(fetch,m) for m in models]):
  receipts.append(future.result());save()
  if len(receipts)%20==0:print(json.dumps({'completed':len(receipts),'total':len(models),'successful':sum(r['success'] for r in receipts)}),flush=True)
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':len(models),'all_requests_pass':all(r['success'] for r in receipts),'empty_models':sum(r.get('records')==0 for r in receipts),'association_tests_run':False}
(p/'metadata/TARGET_RNA_COMPLETE_269.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True)
