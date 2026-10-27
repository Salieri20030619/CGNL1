
from pathlib import Path
import subprocess
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,gzip,datetime,hashlib,select
import numpy as np,pandas as pd
from scipy import io
from threadpoolctl import threadpool_limits
p=Path.cwd();m=p/'mechanism_investigation_v2'
launch=json.loads((m/'metadata/GSE286057_RECOVERY_LAUNCH_209.json').read_text(encoding='utf-8'))
receipt=Path(launch['receipt'])
if not receipt.exists():
 try:
  subprocess.run(['/usr/bin/python3',str(p/'scripts/cgnl_wait_pidfd_system.py'),str(launch['wrapper_pid'])],check=True)
 except ProcessLookupError:pass
assert receipt.exists() and json.loads(receipt.read_text(encoding='utf-8'))['exit_code']==0,'Downloader wrapper did not complete successfully'
assert json.loads((m/'metadata/GSE286057_DOWNLOAD_RECOVERED_COMPLETE_209.json').read_text(encoding='utf-8'))['all_pass']
spec=json.loads((m/'GSE286057_COMPARATOR_SPECIFICATION_195.json').read_text(encoding='utf-8'));panels=spec['source_panels']
clinical=json.loads((m/'metadata/GSE286057_CLINICAL_QUALIFICATION_197.json').read_text(encoding='utf-8'));clin={r['sample_name']:r for r in clinical['samples']}
dest=m/'processed/GSE286057_comparator_198';dest.mkdir(exist_ok=True)
out=m/'results/GSE286057_comparator_198';out.mkdir(exist_ok=True)
needed=sorted(set(sum(panels.values(),[])+['CGNL1','CGN','CAMSAP3','FN1','TNC','ACTA2','TAGLN','KRT8','KRT18','KRT19','PAX8','ITLN1']))
source=[];qc=[];contexts=[];inputs=[]
for sample in spec['samples']:
 title=sample['title'];acc=sample['accession'];root=m/'raw/GSE286057_comparator_195';prefix=acc+'_'+title
 files={k:root/(prefix+'_'+k) for k in ['features.tsv.gz','barcodes.tsv.gz','matrix.mtx.gz']}
 ft=pd.read_csv(files['features.tsv.gz'],sep='	',header=None,dtype=str);bc=pd.read_csv(files['barcodes.tsv.gz'],sep='	',header=None,dtype=str)[0];genes=ft[1].to_numpy();ids=ft[0].to_numpy()
 assert sum(ids=='ENSG00000128849')==1 and genes[np.flatnonzero(ids=='ENSG00000128849')[0]]=='CGNL1'
 with gzip.open(files['matrix.mtx.gz'],'rb') as f:
  with threadpool_limits(limits=2):X=io.mmread(f).tocsr()
 X.sum_duplicates();X.eliminate_zeros()
 assert X.shape==(len(ft),len(bc)) and bc.is_unique and (X.data>=0).all() and np.equal(X.data,np.floor(X.data)).all()
 umi=np.asarray(X.sum(0)).ravel();ng=np.diff(X.tocsc().indptr)
 mito=np.asarray(X[np.char.startswith(genes.astype(str),'MT-')].sum(0)).ravel()/np.maximum(umi,1)
 c={g:np.asarray(X[genes==g].sum(0)).ravel() for g in needed}
 def det(gs):return sum((c[g]>0).astype(int) for g in gs)
 immune=det(panels['immune'])>=2
 epi_flag=(c['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)
 meso=(det(panels['mesothelial_specific'])>=1)&(det(panels['mesothelial_support'])>=1)&~immune
 epi=epi_flag&(det(panels['mesothelial_specific'])==0)&~immune
 fib=(det(panels['fibroblast'])>=3)&~immune&~meso&~epi
 end=(det(panels['endothelial'])>=2)&~immune&~meso&~epi
 mural=(det(panels['mural'])>=2)&~immune&~meso&~epi
 overlap=fib.astype(int)+end.astype(int)+mural.astype(int)>1
 groups={'immune_like':immune,'strict_mesothelial':meso,'conservative_epithelial':epi,'fibroblast_like':fib&~overlap,'endothelial_like':end&~overlap,'mural_like':mural&~overlap,'overlapping_stromal':overlap}
 groups['unclassified']=~np.logical_or.reduce(list(groups.values()));assert np.all(np.stack(list(groups.values())).sum(0)==1)
 groups.update(mesothelial_both_specific=meso&(det(panels['mesothelial_specific'])==2),mesothelial_twoUMI=meso&((c['UPK3B']+c['LRRN4'])>=2)&(sum(c[g] for g in panels['mesothelial_support'])>=2),mesothelial_no_epithelial_flag=meso&~epi_flag)
 meta={k:sample[k] for k in ['accession','title','donor','site','cohort']};meta.update(diagnosis=clin[title]['diagnosis'],biopsy_site=clin[title]['biopsy_site'])
 per=pd.DataFrame(dict(barcode=bc,UMI=umi,n_genes=ng,mito_fraction=mito,epithelial_flag=epi_flag,**c,**groups))
 for scheme in ['primary','sensitivity']:
  q=spec['qc_'+scheme];keep=(umi>=q['UMI_min'])&(ng>=q['genes_min'])&(ng<=q['genes_max'])&(mito<=q['mitochondrial_fraction_max']);per['qc_'+scheme]=keep
  qc.append(dict(meta,scheme=scheme,raw_cells=len(bc),retained_cells=int(keep.sum()),target_UMI=int(c['CGNL1'][keep].sum())))
  for group,mask0 in groups.items():
   mask=mask0&keep;n=int(mask.sum());total=int(umi[mask].sum());target=c['CGNL1'][mask]
   source.append(dict(meta,scheme=scheme,group=group,n_cells=n,target_UMI=int(target.sum()),total_UMI=total,target_detection=float((target>0).mean()) if n else np.nan,log2CPM1=float(np.log2(1+1e6*target.sum()/total)) if total else np.nan,cell_fraction=n/int(keep.sum()) if keep.any() else np.nan,epithelial_flag_cells=int((epi_flag&mask).sum()),epithelial_flag_fraction=float(epi_flag[mask].mean()) if n else np.nan,qualified_cell_count=n>=20))
   if group.startswith('mesothelial') or group=='strict_mesothelial':
    for gene in needed:contexts.append(dict(meta,scheme=scheme,group=group,gene=gene,UMI=int(c[gene][mask].sum()),total_UMI=total,log2CPM1=float(np.log2(1+1e6*c[gene][mask].sum()/total)) if total else np.nan))
 per.to_csv(dest/(prefix+'_cells.tsv.gz'),sep='	',index=False,encoding='utf-8')
 inputs.append(dict(sample=title,files={k:{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for k,f in files.items()}))
 print('SOURCE',title,'QC',qc[-2]['retained_cells'],flush=True)
pd.DataFrame(source).to_csv(out/'source_profiles.tsv',sep='	',index=False,encoding='utf-8')
pd.DataFrame(qc).to_csv(out/'QC.tsv',sep='	',index=False,encoding='utf-8')
pd.DataFrame(contexts).to_csv(out/'context_genes.tsv',sep='	',index=False,encoding='utf-8')
(out/'INPUTS.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
(out/'SOURCE_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':36,'donors':15,'all_pass':True,'statistics_pending':True,'independent_QA_pending':True,'identity':'marker-defined; genomic identity not established'},indent=2),encoding='utf-8')
