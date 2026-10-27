from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,gzip,tarfile,datetime
import pandas as pd,numpy as np
from scipy import io
from threadpoolctl import threadpool_limits
O=Path.cwd()/'mechanism_investigation_v2';ready=O/'metadata/GSE308231_DOWNLOAD.json'
def live():
 for d in Path('/proc').iterdir():
  if d.name.isdigit():
   try:
    if any(a.endswith(b'/cgnl_mech04_download_tissue.py') for a in (d/'cmdline').read_bytes().split(b'\0')):return True
   except (FileNotFoundError,PermissionError):pass
 return False
while not ready.exists() and live():time.sleep(900)
assert ready.exists(),'tissue downloader exited without receipt';receipt=json.loads(ready.read_text());assert receipt['success'],receipt
spec=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text());panels=spec['source_panels'];rows=[];qc=[];needed=set(sum(panels.values(),[])+['CGNL1','KRT8','KRT18','KRT19']);out=O/'processed/GSE308231';out.mkdir(parents=True,exist_ok=True)
with tarfile.open(O/'raw/GSE308231_RAW.tar') as tf:
 for sample in spec['samples']:
  acc=sample['accession'];prefix=next(n[:-len('_matrix.mtx.gz')] for n in tf.getnames() if n.startswith(acc) and n.endswith('_matrix.mtx.gz'))
  def rd(suffix):return gzip.GzipFile(fileobj=tf.extractfile(prefix+suffix))
  ft=pd.read_csv(rd('_features.tsv.gz'),sep='\t',header=None);bc=pd.read_csv(rd('_barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str);genes=ft[1].astype(str).to_numpy();assert (genes=='CGNL1').sum()==1
  with threadpool_limits(limits=2):X=io.mmread(rd('_matrix.mtx.gz')).tocsr()
  assert X.shape==(len(genes),len(bc)) and bc.is_unique and (X.data>=0).all() and np.equal(X.data,np.floor(X.data)).all()
  umi=np.asarray(X.sum(0)).ravel();ng=np.asarray((X>0).sum(0)).ravel();mito=np.asarray(X[np.char.startswith(genes.astype(str),'MT-')].sum(0)).ravel()/np.maximum(umi,1);counts={g:np.asarray(X[genes==g].sum(0)).ravel() for g in needed if g in genes}
  def det(gs):return sum((counts[g]>0).astype(int) for g in gs if g in counts)
  immune=det(panels['immune'])>=2;meso=(det(panels['mesothelial_specific'])>=1)&(det(panels['mesothelial_support'])>=1)&~immune;epi=(counts['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)&(det(panels['mesothelial_specific'])==0)&~immune
  fib=(det(panels['fibroblast'])>=3)&~immune&~meso&~epi;end=(det(panels['endothelial'])>=2)&~immune&~meso&~epi;mural=(det(panels['mural'])>=2)&~immune&~meso&~epi
  overlap=(fib.astype(int)+end.astype(int)+mural.astype(int))>1
  groups={'immune_like':immune,'strict_mesothelial':meso,'conservative_epithelial':epi,'fibroblast_like':fib&~overlap,'endothelial_like':end&~overlap,'mural_like':mural&~overlap,'overlapping_stromal':overlap}
  groups['unclassified']=~np.logical_or.reduce(list(groups.values()));assert np.all(np.stack(list(groups.values())).sum(0)==1)
  typ='peritoneal_metastasis' if any('peritoneal metastasis' in z for z in sample['characteristics_ch1']) else 'primary'
  for scheme in ['primary','sensitivity']:
   q=spec['qc_'+scheme];keep=(umi>=q['UMI_min'])&(ng>=q['genes_min'])&(mito<=q['mitochondrial_fraction_max'])
   if 'genes_max' in q:keep&=ng<=q['genes_max']
   qc.append({'sample':acc,'type':typ,'scheme':scheme,'raw_cells':len(bc),'retained_cells':int(keep.sum()),'target_UMI':float(counts['CGNL1'][keep].sum())})
   for group,mask in groups.items():
    mask=mask&keep;n=int(mask.sum());total=float(umi[mask].sum());target=counts['CGNL1'][mask];rows.append({'sample':acc,'type':typ,'scheme':scheme,'group':group,'n_cells':n,'target_UMI':float(target.sum()),'total_UMI':total,'target_detected_cells':int((target>0).sum()),'target_detection':float((target>0).mean()) if n else np.nan,'CGNL1_ln1pCP10k':float(np.log1p(target.sum()/total*1e4)) if total else np.nan})
  print('SOURCE',acc,typ,qc[-2:],flush=True)
pd.DataFrame(rows).to_csv(O/'results/GSE308231_source_profiles.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(qc).to_csv(O/'results/GSE308231_QC.tsv',sep='\t',index=False,encoding='utf-8')
(O/'results/GSE308231_SOURCE_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(spec['samples']),'status':'descriptive_marker_source','no_inferential_statistics':True,'identity_not_malignancy_confirmation':True,'doublet_filtering_not_reproduced':True},indent=2),encoding='utf-8');print('Tissue source qualification completed',flush=True)
