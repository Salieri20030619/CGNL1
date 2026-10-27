import pathlib,tarfile,gzip,json,numpy as np,pandas as pd
from scipy import io
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'replacement_candidates';gs=json.loads((O/'cell_source_screen_plan.json').read_text())['genes'];rows=[];q=[]
with tarfile.open(V/'02_raw/GSE228598_RAW.tar') as tar:
 for prefix in sorted(n[:-len('_matrix.mtx.gz')] for n in tar.getnames() if n.endswith('_matrix.mtx.gz')):
  def f(n):return gzip.GzipFile(fileobj=tar.extractfile(n))
  ft=pd.read_csv(f(prefix+'_features.tsv.gz'),sep='\t',header=None);sy=ft.iloc[:,1].astype(str).values;bc=pd.read_csv(f(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0]
  with threadpool_limits(limits=2):X=io.mmread(f(prefix+'_matrix.mtx.gz')).tocsr()
  m=pd.read_csv(V/'03_processed'/(prefix+'_gating_refined.tsv.gz'),sep='\t');ix=pd.Index(bc).get_indexer(m.barcode);assert (ix>=0).all();X=X[:,ix];gi=np.flatnonzero(np.isin(sy,gs));Y=X[gi];gg=sy[gi];assert len(set(gg))==len(gg);umi=np.asarray(X.sum(axis=0)).ravel()
  for label,ix in m.groupby('label').groups.items():
   ix=np.array(list(ix));cnt=np.asarray(Y[:,ix].sum(axis=1)).ravel();det=np.asarray((Y[:,ix]>0).sum(axis=1)).ravel();lib=umi[ix].sum()
   for k,g in enumerate(gg):rows.append({'patient':m.patient.iloc[0],'sample':prefix,'sample_type':m.sample_type.iloc[0],'celltype':label,'gene':g,'n_cells':len(ix),'n_detected':int(det[k]),'fraction_detected':det[k]/len(ix),'CPM':cnt[k]/max(lib,1)*1e6})
  q.append({'sample':prefix,'cells':len(m),'genes_found':len(gg)});print(prefix,len(m),flush=True)
a=pd.DataFrame(rows);a['log2CPM']=np.log2(a.CPM+1);a.to_csv(O/'independent_candidate_compartments.tsv.gz',sep='\t',index=False)
e=a[a.n_cells>=20];z=e.groupby(['gene','celltype']).agg(n_patients=('patient','nunique'),detected_patients=('n_detected',lambda v:int((v>0).sum())),median_detection=('fraction_detected','median'),median_log2CPM=('log2CPM','median')).reset_index();z.to_csv(O/'independent_candidate_summary.tsv',sep='\t',index=False)
(O/'INDEPENDENT_COMPLETE.json').write_text(json.dumps({'exit_code':0,'patients':a.patient.nunique(),'samples_read':len(q),'planned_genes':len(gs),'scope':'descriptive localization, NOT independent metastatic state validation; coarse marker-defined compartments'},indent=2));print(z[z.gene.isin(['RAP1GAP2','DSEL','PTPRQ'])].to_string(index=False))
