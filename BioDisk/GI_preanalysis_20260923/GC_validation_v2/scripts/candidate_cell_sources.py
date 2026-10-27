import pathlib,json,numpy as np,pandas as pd
from scipy import io
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_discovery'
s=pd.read_csv(O/'screen_pass_all.tsv',sep='\t');genes=s.loc[(s.n_supported_models>=3)&~s.paused_FABP4,'symbol'].dropna().tolist();(O/'cell_source_screen_plan.json').write_text(json.dumps({'genes':genes,'selection':'all first-screen genes supported in >=3 of 4 sensitivity models, excluding paused FABP4','min_cells':20,'inference':'patient-level descriptive; no cell-level tests'},indent=2))
meta=pd.read_csv(V/'03_processed/GSE234129_annotation_qc.tsv.gz',sep='\t');bar=pd.read_csv(B/'02_raw/GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0];ft=pd.read_csv(B/'02_raw/GSE234129_features.tsv.gz',sep='\t',header=None);sy=ft.iloc[:,1 if ft.shape[1]>1 else 0].astype(str).values
with threadpool_limits(limits=2):X=io.mmread(B/'02_raw/GSE234129_count_matrix.mtx.gz').tocsr()
ix=pd.Index(bar).get_indexer(meta.cell_barcodes);assert (ix>=0).all();X=X[:,ix];sel=np.flatnonzero(np.isin(sy,genes));Y=X[sel];ss=sy[sel];assert len(ss)==len(set(ss));libraries=np.asarray(X.sum(axis=0)).ravel();rows=[]
for (patient,sample,bt),ix in meta.groupby(['patient','sample','broad_type']).groups.items():
 ix=np.asarray(list(ix));counts=np.asarray(Y[:,ix].sum(axis=1)).ravel();det=np.asarray((Y[:,ix]>0).sum(axis=1)).ravel();lib=libraries[ix].sum()
 for k,g in enumerate(ss):rows.append({'patient':patient,'sample':sample,'site':sample.split('-')[-1],'celltype':bt,'gene':g,'n_cells':len(ix),'CPM':counts[k]/max(lib,1)*1e6,'fraction_detected':det[k]/len(ix)})
a=pd.DataFrame(rows);a['log2CPM']=np.log2(a.CPM+1);a.to_csv(O/'candidate_celltype_sample.tsv',sep='\t',index=False)
e=a[a.n_cells>=20];p=e.groupby(['patient','celltype','gene'])[['log2CPM','fraction_detected']].mean().reset_index();q=p.groupby(['gene','celltype']).agg(n_patients=('patient','nunique'),median_patient_log2CPM=('log2CPM','median'),median_patient_detection=('fraction_detected','median')).reset_index();q.to_csv(O/'candidate_cell_source.tsv',sep='\t',index=False)
out=[];paired=[]
for (g,bt),z in e[e.site.isin(['Ca','As'])].groupby(['gene','celltype']):
 w=z.pivot(index='patient',columns='site',values='log2CPM')
 if not {'Ca','As'}<=set(w.columns):continue
 w=w.dropna(subset=['Ca','As']);d=w.As-w.Ca
 if len(d):
  out.append({'gene':g,'celltype':bt,'n_pairs':len(d),'median_delta':d.median(),'positive':sum(d>0),'negative':sum(d<0),'zero':sum(d==0)})
  for pt in d.index:paired.append({'gene':g,'celltype':bt,'patient':pt,'primary_log2CPM':w.loc[pt,'Ca'],'ascites_log2CPM':w.loc[pt,'As'],'delta':d[pt]})
c=pd.DataFrame(out);c.to_csv(O/'candidate_celltype_paired_summary.tsv',sep='\t',index=False);pd.DataFrame(paired).to_csv(O/'candidate_celltype_paired_effects.tsv',sep='\t',index=False)
focus=['COPZ2','TMEM91','CDH23','CPTP','RAP1GAP2','OPLAH','IGSF10','CFI','ACSS3']
print(q[q.gene.isin(focus)&(q.median_patient_detection>=.03)].to_string(index=False));print(c[c.gene.isin(focus)].to_string(index=False));(O/'CELL_SOURCE_COMPLETE.json').write_text(json.dumps({'exit_code':0,'genes_planned':len(genes),'genes_found':len(ss),'cells':len(meta)}))
