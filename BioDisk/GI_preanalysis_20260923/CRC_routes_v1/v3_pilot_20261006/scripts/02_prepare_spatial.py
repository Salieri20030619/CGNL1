from pathlib import Path
import os,json,time,hashlib,datetime
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ.update(NO_PROXY='*',no_proxy='*',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8')
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import anndata as ad
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006')
B=ROOT.parent;BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
plan=json.loads((ROOT/'plan_v3.json').read_text());panels=plan['panels']
def read10x(p):
 with h5py.File(p) as h:
  g=h['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr()
  names=[v.decode() for v in g['features/name'][:]];barcodes=[v.decode() for v in g['barcodes'][:]]
 if len(set(names))<len(names):
  uniq=sorted(set(names));idx={s:i for i,s in enumerate(uniq)}
  X=X@sparse.csr_matrix((np.ones(len(names)),(range(len(names)),[idx[s] for s in names])),shape=(len(names),len(uniq)));names=uniq
 return X.astype(np.float32),np.array(names),barcodes
def scores(X,names):
 umi=np.asarray(X.sum(axis=1)).ravel();norm=sparse.diags(1e4/np.maximum(umi,1))@X;norm=norm.tocsr();norm.data=np.log1p(norm.data)
 out={};cov={}
 for k,gs in panels.items():
  ix=np.flatnonzero(np.isin(names,gs));cov[k]=len(ix)
  out[k]=np.asarray(norm[:,ix].mean(axis=1)).ravel() if len(ix)>=2 else np.full(X.shape[0],np.nan)
 return pd.DataFrame(out),umi,norm,cov
def standard(x):
 x=np.asarray(x,dtype=float);s=x.std()
 return (x-x.mean())/s if s>0 else np.zeros(len(x))
rows=[];geom=[];effects=[];technical=[];coverage=[];allprofs=[];plots=[]
folder=B/'expanded_screen_20260925/spatial_data/GSE217414'
for num,p in enumerate(sorted(folder.glob('*_filtered_feature_bc_matrix.h5'))):
 sample=p.name.split('_')[0];X,names,bc=read10x(p)
 sp=next((folder/(sample+'_unpacked')).glob('*spatial'))
 pos=pd.read_csv(sp/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col']).set_index('barcode')
 assert len(pos.index.unique())==len(pos)
 assert set(bc)<=set(pos.index)
 pos=pos.loc[bc];scale=json.loads((sp/'scalefactors_json.json').read_text())
 xy=pos[['px_col','px_row']].to_numpy(float)*55/scale['spot_diameter_fullres']
 umi=np.asarray(X.sum(axis=1)).ravel();keep=(umi>=500)&(pos.in_tissue.to_numpy()==1)
 X=X[keep];xy=xy[keep];bc=np.array(bc)[keep]
 sc,umi,norm,cov=scores(X,names);target_index=np.flatnonzero(names=='CGNL1');assert len(target_index)==1
 target_raw=X[:,target_index[0]].toarray().ravel();target=norm[:,target_index[0]].toarray().ravel()
 sc['CGNL1']=target;sc['target_raw']=target_raw;sc['UMI']=umi;sc['sample']=sample;sc['patient']=f'GSE217414:{num+1}';sc['barcode']=bc
 tree=cKDTree(xy);dist,_=tree.query(xy,k=2);pitch=float(np.median(dist[:,1]))
 for radius in [100,200,300]:
  pairs=tree.query_pairs(radius,output_type='ndarray');r=np.r_[pairs[:,0],pairs[:,1]];c=np.r_[pairs[:,1],pairs[:,0]]
  A=sparse.csr_matrix((np.ones(len(r)),(r,c)),shape=(len(xy),len(xy)))
  deg=np.asarray(A.sum(axis=1)).ravel();W=sparse.diags(1/np.maximum(deg,1))@A
  mural=np.asarray(W@sc['mural'].values);tcell=np.asarray(W@sc['tcell'].values)
  cols=[mural,sc.endothelial_identity.values,np.log1p(umi),sc.fibroblast.values,sc.epithelial.values,sc.hepatocyte.values]
  design=np.column_stack([np.ones(len(xy))]+[standard(v) for v in cols])
  beta=np.linalg.lstsq(design,standard(target),rcond=None)[0]
  residual=standard(target)-design@beta
  effects.append({'sample':sample,'patient':sc.patient.iloc[0],'radius_um':radius,'spots':len(xy),'beta_neighbor_mural_adjusted':float(beta[1]),'P':None,'FDR':None,'reason_no_P':'technical_mixed_spot_pilot_no_patient_inference','target_tcell_rho':float(spearmanr(target,tcell).statistic),'condition_number':float(np.linalg.cond(design))})
  geom.append({'sample':sample,'radius_um':radius,'median_pitch_um':pitch,'min_neighbors':int(deg.min()),'median_neighbors':float(np.median(deg)),'isolated_spots':int((deg==0).sum()),'cross_section_edges':0})
  if radius==200:
   sc['neighbor_mural']=mural;sc['neighbor_Tcell']=tcell;sc['x_um']=xy[:,0];sc['y_um']=xy[:,1]
   sc.to_csv(ROOT/'results'/f'{sample}_spot_profiles.tsv.gz',sep='\t',index=False)
   if num<2:
    # Geometric contiguous window, target-independent.
    center=np.median(xy,axis=0);sel=np.argsort(np.sum((xy-center)**2,axis=1),kind='stable')[:1500]
    obs=sc.iloc[sel].copy();obs.index=[sample+':'+s for s in obs.barcode]
    a=ad.AnnData(X=X[sel].copy(),obs=obs,var=pd.DataFrame(index=names))
    a.obsm['spatial']=xy[sel];a.obsp['spatial_connectivities']=A[sel][:,sel].tocsr()
    a.layers['counts']=a.X.copy();a.uns['measurement']='mixed Visium spot counts'
    a.write_h5ad(BIG/'inputs'/f'{sample}_technical.h5ad',compression='gzip')
    technical.append({'sample':sample,'path':str(BIG/'inputs'/f'{sample}_technical.h5ad'),'n_spots':len(sel),'selection':'closest to geometric median; no target selection'})
 rows.append({'dataset':'GSE217414','sample':sample,'patient':sc.patient.iloc[0],'n_spots_raw':len(keep),'n_spots_QC':len(xy),'features':len(names),'CGNL1_measured':True,'CGNL1_counts':float(target_raw.sum()),'CGNL1_positive_spots':int((target_raw>0).sum()),'integer_counts':bool(np.allclose(X.data,np.round(X.data))),'pathology_labels_available':False,'cell_resolution':False})
 coverage.append({'sample':sample,**cov})
 plots.append((sample,xy,target,sc.neighbor_mural.values,sc.neighbor_Tcell.values))
pd.DataFrame(rows).to_csv(ROOT/'results/input_qualification.tsv',sep='\t',index=False)
pd.DataFrame(geom).to_csv(ROOT/'results/spatial_geometry.tsv',sep='\t',index=False)
pd.DataFrame(effects).to_csv(ROOT/'results/mixed_spot_effects_descriptive.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(ROOT/'results/panel_coverage.tsv',sep='\t',index=False)
(ROOT/'results/TECHNICAL_INPUTS.json').write_text(json.dumps(technical,indent=2))
fig,axs=plt.subplots(4,3,figsize=(11,13))
for i,(sam,xy,y,m,t) in enumerate(plots):
 for j,(v,label) in enumerate([(y,'CGNL1 log1p CP10k'),(m,'Mural-neighbor score'),(t,'T-cell-neighbor score')]):
  ax=axs[i,j];z=ax.scatter(xy[:,0],-xy[:,1],c=v,s=3,cmap='viridis');ax.set_title(sam+' | '+label,fontsize=9);ax.set_aspect('equal');ax.axis('off');fig.colorbar(z,ax=ax,shrink=.65)
fig.suptitle('Real-data technical pilot: mixed spots; four patient-labelled sections',fontsize=12)
fig.tight_layout();fig.savefig(ROOT/'results/spatial_pilot_maps.png',dpi=160);plt.close(fig)
# Reuse frozen accepted-EC count caches for new panel measurement only.
reg=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text())
lib=[]
for r in reg:
 X=sparse.load_npz(r['path']).tocsr();gs=pd.read_csv(r['genes'],sep='\t')
 # Adapt explicit gene table column only; never select by target outcome.
 names=gs.iloc[:,0].astype(str).to_numpy()
 if X.shape[1]!=len(names) and X.shape[0]==len(names):X=X.T.tocsr()
 if X.shape[1]!=len(names):raise ValueError(('gene dimensions',r['sample'],X.shape,gs.shape))
 total=np.asarray(X.sum(axis=0)).ravel();n=X.shape[0]
 profile=np.log1p(total/max(total.sum(),1)*1e4)
 z={'dataset':r['dataset'],'sample':r['sample'],'patient':r['donor'],'n_accepted_EC':n,'total_UMI':float(total.sum()),'eligible_v3_n20':bool(n>=20)}
 for k,panel in {'CGNL1':['CGNL1'],**panels}.items():
  ix=np.flatnonzero(np.isin(names,panel));z[k]=float(profile[ix].mean()) if len(ix)>=(1 if k=='CGNL1' else 2) else None;z[k+'_covered']=len(ix)
 lib.append(z)
libs=pd.DataFrame(lib);libs.to_csv(ROOT/'results/accepted_EC_library_measurements.tsv',sep='\t',index=False)
elig=libs[libs.eligible_v3_n20].copy()
if len(elig):
 elig.groupby(['dataset','patient'],as_index=False)[['CGNL1','notch_response','adhesion','hypoxia','historical_junction']].mean().to_csv(ROOT/'results/accepted_EC_patient_measurements.tsv',sep='\t',index=False)
done={'exit_code':0,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spatial_patients':4,'spatial_spots_QC':sum(r['n_spots_QC'] for r in rows),'technical_objects':len(technical),'accepted_EC_libraries':len(libs),'accepted_EC_patients_n20':int(elig.patient.nunique()),'inference_P_values_run':False,'boundary_analysis':'not_run_no_independent_pathology_labels','interpretation':'measurement_and_descriptive_baseline_only'}
(ROOT/'results/02_COMPLETE.json').write_text(json.dumps(done,indent=2));print(json.dumps(done))
