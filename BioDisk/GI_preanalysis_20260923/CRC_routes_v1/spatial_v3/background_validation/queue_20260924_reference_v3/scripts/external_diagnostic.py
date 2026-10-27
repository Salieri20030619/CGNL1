from common import *
from mixlib import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
D=Q/'V04_DIAGNOSTIC';D.mkdir(exist_ok=True);rows=[]
for name in ['patient4','patient5']:
 A=Q/'V04_INPUT_QC';X=sparse.load_npz(A/(name+'_counts.npz'));z=np.load(A/(name+'_features.npz'));sy=z['genes'];bc=z['barcodes'];pos=pd.read_csv(A/(name+'_positions.tsv'),sep='\t').set_index('barcode').reindex(bc);lib=np.asarray(X.sum(axis=0)).ravel();good=(lib>=500)&(X.getnnz(axis=0)>=200);types,frac,pred,diag=fit_reference(X,sy,lib);ct=list(types);epi=frac[:,[j for j,t in enumerate(ct) if t in ['epithelial_CNV_supported','epithelial_unresolved']]].sum(axis=1);hep=frac[:,ct.index('hepatocyte')];cnv=frac[:,ct.index('epithelial_CNV_supported')];fib=frac[:,ct.index('fibroblast')]
 row={'sample':name,'QC_spots':int(good.sum()),'hep70':int((good&(hep>=.7)).sum()),'CNV60':int((good&(cnv>=.6)).sum()),'epiTotal60':int((good&(epi>=.6)).sum()),'original_tumor_like':int((good&(cnv>=.6)&(hep<.1)&(fib<.2)).sum()),'broader_epithelial_like':int((good&(epi>=.6)&(hep<.1)&(fib<.2)).sum())};rows.append(row)
 v=pos.reset_index();v['QC']=good
 for j,t in enumerate(types):v['RNAweight_'+t]=frac[:,j]
 v.to_csv(D/(name+'_weights.tsv.gz'),sep='\t',index=False)
 scal=json.loads((A/(name+'_scalefactors.json')).read_text());im=np.asarray(Image.open(A/(name+'_tissue_hires_image.png')));px=pos.pxl_col_in_fullres.to_numpy()*scal['tissue_hires_scalef'];py=pos.pxl_row_in_fullres.to_numpy()*scal['tissue_hires_scalef'];cg=X[np.flatnonzero(sy=='CGNL1')[0]].toarray().ravel()/np.maximum(lib,1)*1e6
 fig,axs=plt.subplots(2,3,figsize=(13,9))
 for ax,(title,value) in zip(axs.ravel(),[('H&E',None),('Hepatic RNA weight',hep),('Total epithelial RNA weight',epi),('CNV-supported epithelial RNA weight',cnv),('Fibroblast RNA weight',fib),('CGNL1 log2(CPM+1)',np.log2(1+cg))]):
  ax.imshow(im)
  if value is not None:
   sc=ax.scatter(px[good],py[good],c=value[good],s=12,alpha=.8,cmap='viridis',vmin=0,vmax=1 if 'weight' in title else None);fig.colorbar(sc,ax=ax,shrink=.65)
  ax.set_title(title);ax.axis('off')
 fig.suptitle(name+' | diagnostic overlays; broad epithelial != confirmed malignant');fig.tight_layout();fig.savefig(D/(name+'_diagnostic.png'),dpi=150);plt.close(fig)
pd.DataFrame(rows).to_csv(D/'region_coverage.tsv',sep='\t',index=False)
dump(D/'RESULT.json',{'status':'diagnostic_completed_pending_image_review','coverage':rows,'scope':'No threshold changed; total epithelial sensitivity must be reported separately from CNV-supported primary'})
