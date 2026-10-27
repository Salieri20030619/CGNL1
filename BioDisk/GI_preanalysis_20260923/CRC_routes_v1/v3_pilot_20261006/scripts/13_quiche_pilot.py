from pathlib import Path
import os,json,gzip,io,importlib.util,time,hashlib,datetime
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NO_PROXY='*',no_proxy='*')
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');B=R.parent;D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
src=R/'vendor/quiche_graph.py';spec=importlib.util.spec_from_file_location('quiche_native_graph',src);qu=importlib.util.module_from_spec(spec);spec.loader.exec_module(qu)
cfg={'method':'QUICHE official unmodified graph.py spatial_niches_khop','scope':'native_spatial_module_only_not_full_differential_enrichment','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'python_note':'standalone source module run in Python3.11; full-package metadata declares <3.11','selection':'first sample GSM9815454; up to6000 QC cells nearest geometric median','QC':'transcript_counts>=50 and nucleus_count==1','labels':'provisional independent marker-panel assignment; unreviewed','label_rule':'at least2 expressed markers, top score>=0.5 and top-second margin>=0.1; otherwise unclassified','scales_um':[100,200,300],'n_neighbors':30,'seed':20261006,'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inference_P_values':False}
(R/'results/QUICHE_PARAMETERS.json').write_text(json.dumps(cfg,indent=2))
with gzip.open(D/'downloads/GSM9815454_cells.parquet.gz','rb') as h:cells=pd.read_parquet(io.BytesIO(h.read())).set_index('cell_id')
with h5py.File(B/'spatial_v3/raw/GSM9815454_tumor1_M242T_cell_feature_matrix.h5') as f:
 g=f['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr().astype(np.float32);names=np.array([x.decode() for x in g['features/name'][:]]);ids=[x.decode() for x in g['barcodes'][:]]
assert cells.index.is_unique and set(ids)<=set(cells.index)
cells=cells.loc[ids];mask=(cells.transcript_counts>=50)&(cells.nucleus_count==1)
X=X[mask.values];cells=cells[mask].copy();xy=cells[['x_centroid','y_centroid']].to_numpy()
order=np.argsort(np.sum((xy-np.median(xy,axis=0))**2,axis=1),kind='stable')[:6000];X=X[order];cells=cells.iloc[order].copy();xy=xy[order]
panels=json.loads((R/'plan_v3.json').read_text())['panels'];keys=['endothelial_identity','mural','fibroblast','epithelial','hepatocyte','tcell','myeloid'];scores=[];hits=[];coverage={}
umi=np.asarray(X.sum(axis=1)).ravel();norm=(sparse.diags(1e4/np.maximum(umi,1))@X).tocsr();norm.data=np.log1p(norm.data)
for k in keys:
 ix=np.flatnonzero(np.isin(names,panels[k]));coverage[k]=[str(names[i]) for i in ix]
 scores.append(np.asarray(norm[:,ix].mean(axis=1)).ravel() if len(ix)>=2 else np.zeros(len(cells)))
 hits.append(np.asarray((X[:,ix]>0).sum(axis=1)).ravel())
S=np.array(scores).T;H=np.array(hits).T;best=S.argmax(axis=1);ranked=np.sort(S,axis=1);labels=np.array(keys,dtype=object)[best]
good=(ranked[:,-1]>=.5)&((ranked[:,-1]-ranked[:,-2])>=.1)&(H[np.arange(len(cells)),best]>=2);labels[~good]='unclassified'
obs=cells.copy();obs['candidate_identity']=pd.Categorical(labels);obs['sample']='GSM9815454';obs['patient']='M242'
a=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=names));a.obsm['spatial']=xy
a.write_h5ad(D/'inputs/GSM9815454_quiche_technical.h5ad',compression='gzip')
assert not any(x.upper() in ['CGNL1','JACOP','PCING','ENSG00000128849'] for x in names)
results=[]
for rad in cfg['scales_um']:
 t=time.time();niches,nn=qu.spatial_niches_khop(a,radius=rad,k=10,khop=3,labels_key='candidate_identity',fov_key='sample',min_cell_threshold=3,n_jobs=8)
 assert np.isfinite(niches.to_numpy()).all() and np.allclose(niches.sum(axis=1),1)
 # independent numerical check of 20 rows against the documented kNN/radius definition
 dist,ix=NearestNeighbors(n_neighbors=30,algorithm='kd_tree').fit(xy).kneighbors()
 for i in range(min(20,len(cells))):
  keep=ix[i][dist[i]<rad]
  if len(keep)<3:continue
  expected=pd.Series(labels[keep]).value_counts(normalize=True).reindex(niches.columns,fill_value=0)
  assert np.allclose(niches.loc[a.obs_names[i]].to_numpy(),expected.to_numpy())
 niches.to_csv(R/'results'/f'QUICHE_niches_radius{rad}.tsv.gz',sep='\t')
 results.append({'radius_um':rad,'cells_retained':len(niches),'identity_classes':niches.shape[1],'elapsed_s':time.time()-t,'row_sums_valid':True,'independent_check_20_cells':True})
 print('SCALE_COMPLETE',results[-1],flush=True)
pd.DataFrame({'provisional_label':labels}).value_counts().rename('cells').reset_index().to_csv(R/'results/QUICHE_label_counts.tsv',sep='\t',index=False)
(R/'results/QUICHE_PANEL_COVERAGE.json').write_text(json.dumps(coverage,indent=2))
out={'exit_code':0,'status':'native_partial_run','module':'spatial_niches_khop','runs':results,'n_patients':1,'CGNL1_status':'not_in_feature_panel_not_zero','labels_reviewed':False,'differential_abundance_run':False,'reason':'one technical sample and provisional labels; not biological CGNL1 validation','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'results/QUICHE_COMPLETE.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
