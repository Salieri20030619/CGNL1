import pathlib,sys,json,os,time,subprocess,hashlib,concurrent.futures,threading,urllib.request,gzip,zipfile,traceback
import numpy as np,pandas as pd
from scipy import io as sio
from threadpoolctl import threadpool_limits
R=pathlib.Path(sys.argv[1]);C=json.loads((R/'CONFIG.json').read_text());D=R/'DATA_PREP';D.mkdir(exist_ok=True);files=json.loads((R/'DOWNLOADS.json').read_text());samples=json.loads((R/'SAMPLES.json').read_text());cfg=C['data_stage'];deadline=time.monotonic()+cfg['deadline_hours']*3600;lock=threading.Lock();committed=0;reserved=0;results=[]
def save(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2));t.replace(p)
def fetch(f):
 global committed,reserved
 p=R/f['relative_path'];p.parent.mkdir(parents=True,exist_ok=True);meta=p.with_suffix(p.suffix+'.complete.json');rec=dict(f);reservation=0;reserved_added=False
 try:
  if time.monotonic()>=deadline:raise TimeoutError('run deadline reached')
  size=f.get('expected_bytes',0)
  if not size:
   try:
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':cfg['proxy']}));h=op.open(urllib.request.Request(f['url'],method='HEAD'),timeout=20);size=int(h.headers.get('Content-Length',0))
   except Exception:pass
  reservation=size or 1024**3
  with lock:
   if committed+reserved+reservation>cfg['total_byte_budget']:raise RuntimeError('run byte budget exceeded')
   reserved+=reservation;reserved_added=True
  if not p.exists():
   remaining=max(30,int((deadline-time.monotonic())/(cfg['retry_count']+1)));timeout=min(cfg['per_file_timeout_seconds'],remaining)
   cmd=['curl','-x',cfg['proxy'],'-L','--fail','--retry',str(cfg['retry_count']),'--retry-all-errors','--retry-delay','2','--connect-timeout','25','--max-time',str(timeout),'--max-filesize',str(1024**3),'-C','-','-o',str(p)+'.part',f['url']]
   with p.with_suffix(p.suffix+'.download.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   if r.returncode:raise RuntimeError('curl exit '+str(r.returncode))
   pathlib.Path(str(p)+'.part').replace(p)
  if size:assert p.stat().st_size==size,'Content-Length mismatch'
  if p.name.endswith('.gz'):
   with gzip.open(p,'rb') as z:
    while z.read(1024*1024):pass
  elif p.suffix=='.zip':
   with zipfile.ZipFile(p) as z:assert z.testzip() is None
  with p.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
  rec.update(ok=True,bytes=p.stat().st_size,sha256=h)
  with lock:committed+=p.stat().st_size
 except Exception as e:rec.update(ok=False,error=str(e))
 finally:
  if reserved_added:
   with lock:reserved-=reservation
 save(meta,rec)
 with lock:
  results.append(rec);save(D/'download_progress.json',{'finished':len(results),'expected':len(files),'successful':sum(x['ok'] for x in results),'validated_bytes':committed,'last_file':f['name']})
 print(f['name'],rec['ok'],rec.get('error',''),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=cfg['max_download_threads']) as ex:downloads=list(ex.map(fetch,files))
save(D/'DOWNLOAD_COMPLETE.json',downloads);okpaths={x['relative_path'] for x in downloads if x['ok']};q=[];targets=[]
programs={'hepatocyte':['ALB','ASGR1','APOA1','TTR'],'epithelial':['EPCAM','KRT8','KRT18','KRT19'],'endothelial':['PECAM1','VWF','KDR','EMCN'],'stromal':['COL1A1','COL1A2','DCN','LUM'],'immune':['PTPRC','LST1','CD3D']}
for sample in samples:
 item={k:v for k,v in sample.items() if k!='files'}
 try:
  assert all(x in okpaths for x in sample['files']),'input incomplete'
  get=lambda suffix:R/next(x for x in sample['files'] if x.endswith(suffix))
  ft=pd.read_csv(get('_features.tsv.gz'),sep='\t',header=None);bc=pd.read_csv(get('_barcodes.tsv.gz'),sep='\t',header=None);assert bc[0].is_unique
  with gzip.open(get('_matrix.mtx.gz'),'rb') as f,threadpool_limits(limits=2):X=sio.mmread(f).tocsr()
  X.eliminate_zeros();assert X.shape==(len(ft),len(bc));assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data));sy=ft[1].astype(str).to_numpy();lib=np.asarray(X.sum(axis=0)).ravel();ng=np.asarray(X.getnnz(axis=0)).ravel();mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;good=(lib>=500)&(ng>=200)&(mt<=25)
  v=pd.DataFrame({'barcode':bc[0],'UMI':lib,'n_genes':ng,'MT_pct':mt,'QC_pass':good});coverage={}
  for gene in C['genes']:
   ix=np.flatnonzero(sy==gene);coverage[gene]=len(ix)
   if len(ix)==1:v[gene+'_counts']=X[ix[0]].toarray().ravel()
  for label,markers in programs.items():
   ix=np.flatnonzero(np.isin(sy,markers));coverage[label]=int(len(ix));v[label+'_score']=np.log1p(X[ix].toarray()/np.maximum(lib,1)*1e4).mean(axis=0) if len(ix)>=2 else np.nan
  out=R/'V02_preparation';out.mkdir(exist_ok=True);v.to_csv(out/(sample['sample']+'_target_marker_QC.tsv.gz'),sep='\t',index=False)
  item.update(status='matrix_validated',cells=X.shape[1],QC_cells=int(good.sum()),genes=X.shape[0],coverage=json.dumps(coverage),identity_status='not_annotated; malignancy_unconfirmed')
  for gene in C['genes']:
   if gene+'_counts' in v and good.any():targets.append({'sample':sample['sample'],'title':sample['title'],'dataset':sample['dataset'],'patient_label_provisional':sample['patient_label_provisional'],'gene':gene,'QC_cells':int(good.sum()),'counts':float(v.loc[good,gene+'_counts'].sum()),'detected_pct':float((v.loc[good,gene+'_counts']>0).mean()*100),'scope':'sample-wide QC only; not cell-source or malignant-cell proof'})
 except Exception as e:item.update(status='failed',error=str(e))
 q.append(item);pd.DataFrame(q).to_csv(D/'sample_QC.tsv',sep='\t',index=False);print(sample['sample'],item['status'],flush=True)
pd.DataFrame(targets).to_csv(D/'sample_target_coverage.tsv',sep='\t',index=False)
archives=[]
for f in files:
 if f['dataset']=='NicheDE' and f['relative_path'] in okpaths:
  with zipfile.ZipFile(R/f['relative_path']) as z:archives.extend({'sample':f['sample'],'archive':f['name'],'member':i.filename,'bytes':i.file_size} for i in z.infolist())
pd.DataFrame(archives).to_csv(D/'external_archive_inventory.tsv',sep='\t',index=False)
success=all(x['ok'] for x in downloads) and all(x['status']=='matrix_validated' for x in q)
save(D/'RESULT_STATUS.json',{'status':'inputs_QC_ready_pending_identity_and_clinical_review' if success else 'partial_failure','files_successful':sum(x['ok'] for x in downloads),'files_total':len(files),'samples_QC_successful':sum(x['status']=='matrix_validated' for x in q),'V02_reference_status':'not_complete; cell identity and clinical review required','V04_validation_status':'not_complete; external pathology/source review required','V03':'dependency_not_met','V05':'conditional_not_started'})
sys.exit(0 if success else 1)
