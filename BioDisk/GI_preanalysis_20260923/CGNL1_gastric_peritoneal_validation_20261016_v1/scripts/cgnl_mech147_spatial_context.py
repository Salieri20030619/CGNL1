from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,h5py,numpy as np,pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
O=Path('C:/codex_test/output/CGNL1_spatial_public_20261016');D=O/'extracted_spatial/data';qa=json.loads((O/'SPATIAL_SCHEMA_146.json').read_text());assert qa['matrix_consistency']['plain_counts_equal_h5'];spec=json.loads((O/'PMGC2026_QC_SPECIFICATION_145.json').read_text());panels=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text())['source_panels']
with h5py.File(D/'filtered_feature_bc_matrix.h5') as h:
 g=h['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).tocsr();bar=[b.decode('utf-8') for b in g['barcodes'][:]];genes=np.array([b.decode('utf-8') for b in g['features/name'][:]])
pos=pd.read_csv(D/'spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode').loc[bar];sf=json.loads((D/'spatial/scalefactors_json.json').read_text());umi=np.asarray(X.sum(0)).ravel();ng=np.asarray((X>0).sum(0)).ravel();mt=np.asarray(X[np.char.startswith(genes,'MT-')].sum(0)).ravel()/np.maximum(umi,1)
needed=sorted(set(sum(panels.values(),[])+['CGNL1']));counts={}
for name in needed:
 ix=np.where(genes==name)[0];assert len(ix)==1,(name,len(ix));counts[name]=X[ix[0]].toarray().ravel()
df=pos.copy();df['total_UMI']=umi;df['detected_genes']=ng;df['mitochondrial_fraction']=mt;df['CGNL1_UMI']=counts['CGNL1'];df['CGNL1_ln1pCP10k']=np.log1p(counts['CGNL1']/np.maximum(umi,1)*1e4);df['primary_keep']=(pos.in_tissue==1)&(umi>0);df['sensitivity_keep']=df.primary_keep&(umi>=500)&(ng>=200)&(mt<=.25)
for group,names in panels.items():
 raw=sum(counts[n] for n in names);df[group+'_UMI']=raw;df[group+'_detected_markers']=sum((counts[n]>0).astype(int) for n in names);df[group+'_ln1pCP10k']=np.log1p(raw/np.maximum(umi,1)*1e4)
for name,c in counts.items():df['gene_'+name+'_UMI']=c
df.to_csv(O/'PMGC2026_spatial_spot_values.tsv.gz',sep='\t',encoding='utf-8');summaries=[];corr=[]
for scheme in ['primary','sensitivity']:
 v=df.loc[df[scheme+'_keep']];summaries.append({'scheme':scheme,'spots':len(v),'independent_lesions':1,'CGNL1_total_UMI':int(v.CGNL1_UMI.sum()),'CGNL1_positive_spots':int((v.CGNL1_UMI>0).sum()),'CGNL1_positive_fraction':float((v.CGNL1_UMI>0).mean()),'CGNL1_UMI_max':int(v.CGNL1_UMI.max()),'UMI_median':float(v.total_UMI.median()),'genes_median':float(v.detected_genes.median()),'meso_specific_positive_spots':int((v.mesothelial_specific_UMI>0).sum()),'both_specific_and_support_spots':int(((v.mesothelial_specific_UMI>0)&(v.mesothelial_support_UMI>0)).sum()),'target_in_both_meso_marker_context_spots':int(v.loc[(v.mesothelial_specific_UMI>0)&(v.mesothelial_support_UMI>0),'CGNL1_UMI'].sum())})
 Z=np.c_[np.ones(len(v)),v.total_UMI.rank().to_numpy(),v.detected_genes.rank().to_numpy()];y=v.CGNL1_ln1pCP10k.rank().to_numpy();ry=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0]
 for group in panels:
  x=v[group+'_ln1pCP10k'].rank().to_numpy();rx=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];corr.append({'scheme':scheme,'panel':group,'within_section_rank_correlation':float(np.corrcoef(x,y)[0,1]),'depth_and_genes_adjusted_rank_correlation':float(np.corrcoef(rx,ry)[0,1]),'biological_n':1,'inferential_P':None})
pd.DataFrame(corr).to_csv(O/'PMGC2026_spatial_context_descriptive.tsv',sep='\t',index=False,encoding='utf-8')
im=np.asarray(Image.open(D/'spatial/tissue_hires_image.png'));xp=df.pixel_col*sf['tissue_hires_scalef'];yp=df.pixel_row*sf['tissue_hires_scalef'];views=[('H&E',None),('CGNL1 UMI','CGNL1_UMI'),('CGNL1 ln(1+CP10k)','CGNL1_ln1pCP10k'),('Mesothelial-specific markers','mesothelial_specific_ln1pCP10k'),('Mesothelial-support markers','mesothelial_support_ln1pCP10k'),('Epithelial markers','epithelial_ln1pCP10k'),('Endothelial markers','endothelial_ln1pCP10k'),('Fibroblast markers','fibroblast_ln1pCP10k'),('Mural markers','mural_ln1pCP10k')]
fig,axes=plt.subplots(3,3,figsize=(12,12))
for ax,(title,col) in zip(axes.ravel(),views):
 ax.imshow(im,alpha=1 if col is None else .6);ax.set_title(title,fontsize=11);ax.axis('off')
 if col:
  keep=df.primary_keep&(df[col]>0);s=ax.scatter(xp[keep],yp[keep],c=df.loc[keep,col],s=5,cmap='magma',vmin=0,vmax=max(float(df[col].max()),1e-9),linewidths=0);fig.colorbar(s,ax=ax,fraction=.038,pad=.01)
fig.suptitle('One PMGC lesion | positive spots over full histology | panels have separate scales',fontsize=12);fig.tight_layout(rect=[0,0,1,.97]);fig.savefig(O/'PMGC2026_spatial_context.png',dpi=180);plt.close(fig)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary':summaries,'correlations':corr,'interpretation_limits':['one section; no inferential biological P values','spot mixtures and normalized correlations do not assign CGNL1 to a cell type','UPK3B/LRRN4 and support-marker context not proof of intact mesothelial surface','H&E inspected before target overlay; broad gland-like and stromal regions visible, continuous mesothelial boundary not independently confirmed','independent numeric QA pending'],'source_panels':panels};(O/'SPATIAL_CONTEXT_147.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
