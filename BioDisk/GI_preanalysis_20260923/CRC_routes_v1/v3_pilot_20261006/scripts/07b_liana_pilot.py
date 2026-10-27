from pathlib import Path
import os,json,datetime,time,traceback
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ.update(NO_PROXY='*',no_proxy='*',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
import liana as li
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
# Literature-motivated pairs selected before viewing spatial pair scores.
pairs=pd.DataFrame([('DLL4','NOTCH1'),('DLL4','NOTCH4'),('JAG1','NOTCH1'),('JAG1','NOTCH4'),('VEGFA','KDR'),('PDGFB','PDGFRB'),('CXCL12','CXCR4'),('TGFB1','TGFBR2')],columns=['ligand','receptor'])
cfg={'method':'LIANA native SpatialBivariate','version':li.__version__,'pairs':pairs.to_dict('records'),'local_metric':'cosine','n_perms':None,'nz_prop':0.01,'connectivity':'native LIANA linear kernel radius 200 microns','seed':20261006,'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'P_values':'not_run; technical scores only','cells_vs_spots':'mixed_visium_spots'}
(R/'results/LIANA_PILOT_PARAMETERS.json').write_text(json.dumps(cfg,indent=2))
rows=[]
for rec in json.loads((R/'results/TECHNICAL_INPUTS.json').read_text()):
 start=time.time();a=ad.read_h5ad(rec['path']);sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 li.utils.spatial_neighbors(a,bandwidth=200,kernel='linear',cutoff=0,max_neighbours=100,standardize=True,key_added='liana')
 key=[k for k in a.obsp if k.startswith('liana')][0]
 avail=pairs[pairs.ligand.isin(a.var_names)&pairs.receptor.isin(a.var_names)].copy()
 out=li.mt.bivariate(a,resource=avail,local_name='cosine',global_name=None,connectivity_key=key,n_perms=None,nz_prop=.01,use_raw=False,verbose=False)
 assert isinstance(out,ad.AnnData)
 out.write_h5ad(D/'inputs'/f'{rec["sample"]}_LIANA_scores.h5ad',compression='gzip')
 dense=out.X.toarray() if sparse.issparse(out.X) else out.X
 assert np.isfinite(dense).all()
 stats=out.var.copy();stats['pair']=out.var_names;stats['mean_score']=np.mean(dense,axis=0);stats['sd_score']=np.std(dense,axis=0);stats['sample']=rec['sample']
 stats.to_csv(R/'results'/f'{rec["sample"]}_LIANA_pair_scores.tsv',sep='\t',index=False)
 # Same values, spatial coordinates permuted. Technical response only, not a biological null test.
 rng=np.random.default_rng(20261006);b=a.copy();b.obsm['spatial']=a.obsm['spatial'][rng.permutation(a.n_obs)]
 li.utils.spatial_neighbors(b,bandwidth=200,kernel='linear',cutoff=0,max_neighbours=100,standardize=True,key_added='liana')
 z=li.mt.bivariate(b,resource=avail,local_name='cosine',global_name=None,connectivity_key=key,n_perms=None,nz_prop=.01,use_raw=False,verbose=False)
 arr=z.X.toarray() if sparse.issparse(z.X) else z.X
 assert list(z.var_names)==list(out.var_names)
 mad=float(np.mean(np.abs(dense-arr)));diff=stats[['pair','sample','mean_score']].copy();diff['permuted_mean_score']=arr.mean(axis=0);diff['mean_abs_location_score_change']=np.mean(np.abs(dense-arr),axis=0);diff.to_csv(R/'results'/f'{rec["sample"]}_LIANA_geometry_control.tsv',sep='\t',index=False)
 rows.append({'sample':rec['sample'],'spots':a.n_obs,'tested_pairs':out.n_vars,'finite_scores':True,'spatial_shuffle_mean_abs_change':mad,'elapsed_s':time.time()-start,'status':'native_run_complete','interpretation':'technical_spatial_pair_scores_not_CGNL1_causality'})
pd.DataFrame(rows).to_csv(R/'results/LIANA_RUNS.tsv',sep='\t',index=False)
(R/'results/LIANA_COMPLETE.json').write_text(json.dumps({'exit_code':0,'method':cfg['method'],'rows':rows,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print(json.dumps(rows),flush=True)
