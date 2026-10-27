import pathlib,json,numpy as np,pandas as pd
from scipy import io
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2';P=json.loads((R/'00_docs'/'external_gating_plan.json').read_text())
ft=pd.read_csv(B/'02_raw'/'GSE234129_features.tsv.gz',sep='\t',header=None);genes=ft.iloc[:,1 if ft.shape[1]>1 else 0].astype(str).values;bar=pd.read_csv(B/'02_raw'/'GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str);m=pd.read_csv(B/'04_results'/'GC'/'GSE234129_cell_qc.tsv',sep='\t').set_index('cell_barcodes').loc[bar].reset_index()
with threadpool_limits(limits=2):X=io.mmread(B/'02_raw'/'GSE234129_count_matrix.mtx.gz').tocsr()
k={name:np.asarray((X[np.isin(genes,gs)]>0).sum(axis=0)).ravel() for name,gs in P['lineage'].items()};gate=(k['myeloid_core']>=2)&(k['macrophage_mature']>=2)&(k['T']<2)&(k['B']<2)&(k['epithelial']<2)&m.basic_pass.values;truth=m.celltype.str.startswith('TAM').values&m.basic_pass.values
m['marker_gate']=gate;m['author_macrophage']=truth;rows=[]
for patient,z in m.groupby('patient'):
 tp=int((z.marker_gate&z.author_macrophage).sum());fp=int((z.marker_gate&~z.author_macrophage).sum());fn=int((~z.marker_gate&z.author_macrophage).sum());rows.append({'patient':patient,'TP':tp,'FP':fp,'FN':fn,'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/(tp+fn) if tp+fn else None})
pd.DataFrame(rows).to_csv(R/'04_results'/'independent_scRNA'/'gating_reference_refined_patient_QC.tsv',sep='\t',index=False);pd.crosstab(m.loc[m.basic_pass,'celltype'],m.loc[m.basic_pass,'marker_gate']).to_csv(R/'04_results'/'independent_scRNA'/'gating_reference_refined_confusion.tsv',sep='\t');print(pd.DataFrame(rows).to_string(index=False));tp=sum(z['TP'] for z in rows);fp=sum(z['FP'] for z in rows);fn=sum(z['FN'] for z in rows);print('OVERALL_PRECISION_RECALL',tp/(tp+fp),tp/(tp+fn))
