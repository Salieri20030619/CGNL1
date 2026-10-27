from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,gzip,datetime
import numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'feasibility_audit_20261008_v1';O=R/'results'
old=B/'review_consensus_retest_20261002_v1';reg=json.loads((old/'CACHE_INDEX.json').read_text());ann={ds:pd.read_csv(old/'results'/f'{ds}_accepted_annotation.tsv.gz',sep='\t') for ds in sorted({r['dataset'] for r in reg})}
rows=[];checks=[]
for rec in sorted(reg,key=lambda r:(r['dataset'],r['sample'])):
 print('ATLAS_START',rec['sample'],flush=True)
 stem=rec['source_matrix'].replace('matrix.mtx.gz','');genes=pd.read_csv(stem+'features.tsv.gz',sep='\t',header=None)[1].astype(str).to_numpy();bc=pd.read_csv(stem+'barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
 with gzip.open(rec['source_matrix'],'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 X.eliminate_zeros();assert X.shape==(len(genes),len(bc))
 ec=pd.read_csv(rec['cells'],sep='\t');ie=pd.Index(bc).get_indexer(ec.barcode);Y=sparse.load_npz(rec['path']).tocsr();diff=X[:,ie]-Y;assert (ie>=0).all() and diff.nnz==0
 a=ann[rec['dataset']];a=a[a['sample']==rec['sample']].copy();ii=pd.Index(bc).get_indexer(a.barcode);assert (ii>=0).all();assert not a.barcode.duplicated().any()
 umi=np.asarray(X[:,ii].sum(0)).ravel();assert np.array_equal(umi,a.UMI_QC)
 target=np.asarray(X[np.flatnonzero(genes=='CGNL1')][:,ii].sum(0)).ravel();all_target=float(target.sum())
 for label,idx in a.reset_index(drop=True).groupby('accepted_major_lineage').indices.items():
  v=target[idx];u=umi[idx]
  rows.append(dict(dataset=rec['dataset'],sample=rec['sample'],patient=rec['donor'],lineage=label,cells=len(idx),CGNL1_counts=float(v.sum()),CGNL1_positive_cells=int((v>0).sum()),detection_fraction=float(np.mean(v>0)),pseudobulk_CPM=float(v.sum()/max(u.sum(),1)*1e6),mean_cell_log1p_CP10k=float(np.log1p(v/np.maximum(u,1)*1e4).mean()),fraction_of_sample_CGNL1_counts=float(v.sum()/all_target) if all_target else np.nan,total_UMI=float(u.sum())))
 checks.append(dict(dataset=rec['dataset'],sample=rec['sample'],QC_cells=len(a),original_cache_cells=len(ec),raw_cache_entries_different=diff.nnz,all_lineage_UMI_match=True,CGNL1_feature_rows=int((genes=='CGNL1').sum())))
 pd.DataFrame(rows).to_csv(O/'source_atlas_library_lineage.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'all27_raw_cache_and_atlas_checks.tsv',sep='\t',index=False)
 del X,Y
 print('ATLAS_DONE',rec['sample'],len(a),flush=True)
t=pd.DataFrame(rows)
pt=t.groupby(['dataset','patient','lineage'],as_index=False).agg(libraries=('sample','size'),cells=('cells','sum'),CGNL1_counts=('CGNL1_counts','sum'),positive_cells=('CGNL1_positive_cells','sum'),mean_logCP10k_equal_library=('mean_cell_log1p_CP10k','mean'),CPM_equal_library=('pseudobulk_CPM','mean'))
pt['detection_fraction']=pt.positive_cells/pt.cells;den=pt.groupby(['dataset','patient']).CGNL1_counts.transform('sum');pt['fraction_of_patient_CGNL1_counts']=pt.CGNL1_counts/den.replace(0,np.nan)
pt.to_csv(O/'source_atlas_patient_lineage.tsv',sep='\t',index=False)
summary=pt.groupby(['dataset','lineage'],as_index=False).agg(patients=('patient','nunique'),total_cells=('cells','sum'),total_CGNL1=('CGNL1_counts','sum'),median_patient_detection=('detection_fraction','median'),median_patient_CPM=('CPM_equal_library','median'),median_patient_target_fraction=('fraction_of_patient_CGNL1_counts','median'))
summary.to_csv(O/'source_atlas_summary.tsv',sep='\t',index=False)
tops=pt.sort_values(['dataset','patient','CPM_equal_library'],ascending=[True,True,False]).groupby(['dataset','patient']).first().reset_index();tops.to_csv(O/'source_atlas_highest_CPM_lineage_per_patient.tsv',sep='\t',index=False)
done={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'raw_matrices':len(checks),'QC_cells':sum(c['QC_cells'] for c in checks),'cache_matches_all':all(c['raw_cache_entries_different']==0 for c in checks),'patients':pt[['dataset','patient']].drop_duplicates().shape[0],'interpretation':'descriptive frozen broad lineages; count contribution and expression per cell distinguished; not malignant status confirmation; no patient inferential tests'}
(O/'SOURCE_ATLAS_COMPLETE.json').write_text(json.dumps(done,indent=2));print(json.dumps(done),flush=True)
