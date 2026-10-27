import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import pathlib,json,numpy as np,pandas as pd
from scipy import sparse,io
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2';O=R/'04_results'/'paired_scRNA';O.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(B/'04_results'/'GC'/'GSE234129_cell_qc.tsv',sep='\t');bar=pd.read_csv(B/'02_raw'/'GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str);ft=pd.read_csv(B/'02_raw'/'GSE234129_features.tsv.gz',sep='\t',header=None);genes=ft.iloc[:,1 if ft.shape[1]>1 else 0].astype(str).values
m=m.set_index('cell_barcodes').loc[bar].reset_index();assert list(m.cell_barcodes)==list(bar)
with threadpool_limits(limits=2):X=io.mmread(B/'02_raw'/'GSE234129_count_matrix.mtx.gz').tocsr()
def broad(t):
 if t.startswith('TAM'):return 'Macrophage'
 if 'Mono_' in t:return 'Monocyte'
 if t.startswith('Mesothelial'):return 'Mesothelial'
 if t.startswith('Endothelial'):return 'Endothelial'
 if t.startswith('Fibroblast'):return 'Fibroblast'
 if t.startswith('VSMC'):return 'VSMC'
 if 'DC' in t or t=='pDCs':return 'Dendritic'
 if t=='Mast':return 'Mast'
 return 'Lymphoid'
m['broad_type']=m.celltype.map(broad);m['site']=m['sample'].str.split('-').str[-1]
# No marker-based malignant relabeling; original author annotations retained.
keep=m.basic_pass.astype(bool).values;X=X[:,keep];m=m[keep].reset_index(drop=True)
m.to_csv(R/'03_processed'/'GSE234129_annotation_qc.tsv.gz',sep='\t',index=False)
markers=['FABP4','CD36','LPL','FABP5','SLC27A1','SLC27A2','SLC27A4','LDLR','ADIPOQ','PLIN1','CIDEC','CIDEA','PPARG','TREM2','APOE','APOC1','LIPA','PLIN2','GPNMB','MARCO','MSR1','C1QA','SPP1','DAB2','VSIG4','FOLR2','MRC1','CD163','CD68','LST1','PTPRC','EPCAM','KRT19','COL1A1','DCN','PECAM1','MSLN','CALB2']
idx={g:np.flatnonzero(genes==g) for g in markers};pb=[];meta=[];rec=[]
for (patient,sample,bt),ix in m.groupby(['patient','sample','broad_type']).groups.items():
 ix=np.array(list(ix));y=np.asarray(X[:,ix].sum(axis=1)).ravel();key=sample+'__'+bt;pb.append(y);meta.append({'pb':key,'patient':patient,'sample':sample,'site':sample.split('-')[-1],'celltype':bt,'n_cells':len(ix),'library_size':y.sum(),'eligible':len(ix)>=20})
 for g,gi in idx.items():
  if len(gi):
   vals=np.asarray(X[gi][:,ix].sum(axis=0)).ravel();rec.append({'patient':patient,'sample':sample,'celltype':bt,'gene':g,'n_cells':len(ix),'fraction_detected':float((vals>0).mean()),'pseudobulk_CPM':float(y[gi].sum()/max(y.sum(),1)*1e6)})
pb=np.column_stack(pb);meta=pd.DataFrame(meta);pd.DataFrame(pb,index=genes,columns=meta.pb).to_csv(R/'03_processed'/'GSE234129_pseudobulk_counts.tsv.gz',sep='\t');meta.to_csv(O/'pseudobulk_design.tsv',sep='\t',index=False);rec=pd.DataFrame(rec);rec.to_csv(O/'marker_detection_patient_sample.tsv',sep='\t',index=False)
# Patient-level source summary only, eligibility >=20 cells; within-patient multiple sites averaged.
e=rec[rec.n_cells>=20].copy();e['log2CPM']=np.log2(e.pseudobulk_CPM+1);p=e.groupby(['patient','celltype','gene'])[['log2CPM','fraction_detected']].mean().reset_index();q=p.groupby(['celltype','gene']).agg(median_patient_log2CPM=('log2CPM','median'),median_patient_detection=('fraction_detected','median'),patients=('patient','nunique')).reset_index();q.to_csv(O/'cell_source_patient_summary.tsv',sep='\t',index=False)
# Matched primary/ascites, same cell type. Patient is independent unit.
contr=[]
for bt in ['Macrophage','Monocyte','Lymphoid']:
 z=meta[(meta.celltype==bt)&meta.eligible&meta.site.isin(['Ca','As'])];ct=pd.crosstab(z.patient,z.site);pairs=ct.index[(ct.get('Ca',0)==1)&(ct.get('As',0)==1)]
 for patient in pairs:
  for g in markers:
   za=rec[(rec.patient==patient)&(rec.celltype==bt)&(rec.gene==g)&rec['sample'].str.endswith('-As')];zp=rec[(rec.patient==patient)&(rec.celltype==bt)&(rec.gene==g)&rec['sample'].str.endswith('-Ca')]
   if len(za)==1 and len(zp)==1:contr.append({'patient':patient,'celltype':bt,'gene':g,'ascites_log2CPM':float(np.log2(za.pseudobulk_CPM.iloc[0]+1)),'primary_log2CPM':float(np.log2(zp.pseudobulk_CPM.iloc[0]+1)),'delta_log2CPM':float(np.log2(za.pseudobulk_CPM.iloc[0]+1)-np.log2(zp.pseudobulk_CPM.iloc[0]+1)),'ascites_cells':int(za.n_cells.iloc[0]),'primary_cells':int(zp.n_cells.iloc[0])})
c=pd.DataFrame(contr);c.to_csv(O/'paired_gene_effects.tsv',sep='\t',index=False);cg=c.groupby(['celltype','gene']).agg(n_pairs=('patient','nunique'),median_difference=('delta_log2CPM','median'),positive_pairs=('delta_log2CPM',lambda x:int((x>0).sum())),negative_pairs=('delta_log2CPM',lambda x:int((x<0).sum()))).reset_index();cg.to_csv(O/'paired_gene_summary.tsv',sep='\t',index=False)
summary={'cells_after_qc':len(m),'patients':m.patient.nunique(),'cell_source_units':'patient-level summaries from author labels; >=20 cells per sample/type','paired_macrophage_patients':int(c[c.celltype=='Macrophage'].patient.nunique()),'malignant_labels':0,'statistical_scope':'paired sample fewer than 10; descriptive effects, no cell-level inferential tests'};(O/'summary.json').write_text(json.dumps(summary,indent=2));(O/'COMPLETE.json').write_text(json.dumps({'exit_code':0}));print(json.dumps(summary));print(cg[(cg.celltype=='Macrophage')&cg.gene.isin(['FABP4','CD36','LPL','TREM2','APOE','LIPA','GPNMB'])].to_string(index=False));print(q[q.gene.isin(['FABP4','CD36','LPL'])].to_string(index=False))
