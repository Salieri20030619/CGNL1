from common import *
from mixlib import *
from scipy.spatial import cKDTree
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Q/'V04_BROAD_EP_SENS';D.mkdir(exist_ok=True);gate=json.loads((Q/'V02_SOURCE/GATE.json').read_text())
if not gate['allowed']:dump(D/'RESULT.json',{'status':'skipped_insufficient_reference','reason':gate});sys.exit(0)
qc=[];effects=[]
for sample in ['patient4','patient5']:
 try:
  A=Q/'V04_INPUT_QC';X=sparse.load_npz(A/(sample+'_counts.npz'));z=np.load(A/(sample+'_features.npz'));sy=z['genes'];bc=z['barcodes'];pos=pd.read_csv(A/(sample+'_positions.tsv'),sep='\t').set_index('barcode').reindex(bc);lib=np.asarray(X.sum(axis=0)).ravel();good=(lib>=500)&(np.asarray(X.getnnz(axis=0)).ravel()>=200);types,frac,pred,diag=fit_reference(X,sy,lib);ct=list(types);hep=good&(frac[:,ct.index('hepatocyte')]>=.7);total_ep=frac[:,[j for j,t in enumerate(ct) if t in ['epithelial_CNV_supported','epithelial_unresolved']]].sum(axis=1);tum=good&(total_ep>=.6)&(frac[:,ct.index('hepatocyte')]<.1)&(frac[:,ct.index('fibroblast')]<.2);assert hep.sum()>=20 and tum.sum()>=20,'insufficient_model_dominant_regions';xy=np.c_[50*pos.array_col,50*np.sqrt(3)*pos.array_row];dist=cKDTree(xy[hep]).query(xy)[0];dd,ii=cKDTree(xy).query(xy,k=7,distance_upper_bound=100+1e-4);valid=np.isfinite(dd[:,1:])&(ii[:,1:]<len(bc));nn=tum[np.minimum(ii[:,1:],len(bc)-1)]&valid;strict=tum&(valid.sum(axis=1)==6)&(nn.sum(axis=1)>=5);near=strict&(dist>100+1e-6)&(dist<=300+1e-6);far=strict&(dist>=600-1e-6);q={'sample':sample,'QC_spots':int(good.sum()),'hepatic_dominant':int(hep.sum()),'tumor_like_dominant':int(tum.sum()),'near_spots':int(near.sum()),'far_spots':int(far.sum()),'interpretation':'broad epithelial model regions; malignant identity unresolved; morphology pending',**diag};qc.append(q)
  vv=pos.reset_index();vv['distance_to_hepatic_um']=dist;vv['near']=near;vv['far']=far
  for k,ctname in enumerate(types):vv['RNAweight_'+ctname]=frac[:,k]
  for gene in C['targets']:
   g=np.flatnonzero(sy==gene)
   if len(g)!=1:continue
   counts=X[g[0]].toarray().ravel();vv[gene+'_counts']=counts
   if near.sum()>=20 and far.sum()>=20:
    e=np.log2(1+counts[near].sum()/lib[near].sum()*1e6)-np.log2(1+counts[far].sum()/lib[far].sum()*1e6);effects.append({'sample':sample,'patient':sample,'gene':gene,'mask':'broad_epithelial_model_strict5','effect':float(e),'scope':'n<=2 descriptive; no inferred colon-specific causation'})
  vv.to_csv(D/(sample+'_spatial.tsv.gz'),sep='\t',index=False);fig,axs=plt.subplots(1,2,figsize=(9,4));px=pos[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();color=np.full(len(bc),'#DDD',dtype=object);color[hep]='#6BAED6';color[near]='#D06B45';color[far]='#6750A1';axs[0].scatter(px[:,0],px[:,1],c=color,s=3);axs[0].set_title('RNA-reference region model');cg=vv.get('CGNL1_counts',pd.Series(np.zeros(len(bc)))).to_numpy();sct=axs[1].scatter(px[:,0],px[:,1],c=np.log2(1+cg/np.maximum(lib,1)*1e6),s=3,cmap='viridis');fig.colorbar(sct,ax=axs[1],shrink=.6);axs[1].set_title('CGNL1 log2(CPM+1)')
  for ax in axs:ax.invert_yaxis();ax.set_aspect('equal');ax.axis('off')
  fig.suptitle(sample+' | exploratory; morphology review pending');fig.tight_layout();fig.savefig(D/(sample+'.png'),dpi=160);plt.close(fig)
 except Exception as e:qc.append({'sample':sample,'status':'not_interpretable','reason':str(e)})
pd.DataFrame(qc).to_csv(D/'QC.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(D/'effects.tsv',sep='\t',index=False);dump(D/'RESULT.json',{'status':'exploratory_pending_morphology_review' if effects else 'insufficient_data','effects':effects,'independent_replicates_max':2,'clinical_site':'CRC; colon/rectum unresolved','no_inferential_P_values':True})
