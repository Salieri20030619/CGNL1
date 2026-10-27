from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,csv,numpy as np,pandas as pd,datetime,hashlib
O=Path('C:/codex_test/output/CGNL1_spatial_public_20261016');F=O/'extracted_spatial/data/filtered_feature_bc_matrix';features=list(csv.reader((F/'features.tsv').open(encoding='utf-8'),delimiter='\t'));bar=[r[0] for r in csv.reader((F/'barcodes.tsv').open(encoding='utf-8'))];target=next(i+1 for i,r in enumerate(features) if r[0]=='ENSG00000128849');panels=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text())['source_panels'];selected={i+1:r[1] for i,r in enumerate(features) if r[1] in set(sum(panels.values(),[])+['CGNL1'])};counts={g:np.zeros(len(bar),dtype=np.int64) for g in selected.values()};tot=np.zeros(len(bar),dtype=np.int64);ng=np.zeros(len(bar),dtype=np.int64);mito=np.zeros(len(bar),dtype=np.int64);mi={i+1 for i,r in enumerate(features) if r[1].startswith('MT-')};n=0
with (F/'matrix.mtx').open(encoding='utf-8') as f:
 for line in f:
  if not line.startswith('%'):shape=list(map(int,line.split()));break
 assert shape[:2]==[len(features),len(bar)]
 for line in f:
  i,c,v=map(int,line.split());assert v>0;c-=1;tot[c]+=v;ng[c]+=1;n+=1
  if i in mi:mito[c]+=v
  if i in selected:counts[selected[i]][c]+=v
assert n==shape[2];v=pd.read_csv(O/'PMGC2026_spatial_spot_values.tsv.gz',sep='\t').set_index('barcode').loc[bar];assert np.array_equal(v.total_UMI,tot) and np.array_equal(v.detected_genes,ng);errors=[]
for g,a in counts.items():assert np.array_equal(v['gene_'+g+'_UMI'],a)
for group,gs in panels.items():assert np.array_equal(v[group+'_UMI'],sum(counts[g] for g in gs))
norm=np.log1p(counts['CGNL1']/np.maximum(tot,1)*1e4);err=float(np.max(abs(norm-v.CGNL1_ln1pCP10k)));assert err<1e-12;assert np.max(abs(mito/np.maximum(tot,1)-v.mitochondrial_fraction))<1e-12
summary=json.loads((O/'SPATIAL_CONTEXT_147.json').read_text(encoding='utf-8'))['summary'];qc=(tot>=500)&(ng>=200)&(mito/np.maximum(tot,1)<=.25);assert np.array_equal(qc,v.sensitivity_keep.to_numpy());assert int(counts['CGNL1'].sum())==summary[0]['CGNL1_total_UMI'];assert int((counts['CGNL1']>0).sum())==summary[0]['CGNL1_positive_spots'];assert int(counts['CGNL1'][qc].sum())==summary[1]['CGNL1_total_UMI'];marker_summary={g:{'UMI':int(counts[g].sum()),'positive_spots':int((counts[g]>0).sum())} for g in panels['mesothelial_specific']+panels['mesothelial_support']};out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'method':'independent line-by-line plainMatrixMarket integer accumulation; compared to HDF5-derived spot table','triplets_checked':n,'spots':len(bar),'all_target_and_marker_counts_equal':True,'whole_library_UMI_and_gene_counts_equal':True,'normalization_max_error':err,'QC_selection_equal':True,'mesothelial_marker_coverage':marker_summary,'statistical_scope':'descriptive one lesion; no biological P or FDR calculated','visual_review':'wholeH&E and9panel full-field overlay inspected; no continuous mesothelial layer independently established'};(O/'SPATIAL_NUMERIC_QA_148.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
