from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
import json,hashlib,datetime,importlib.util,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');A=B/'annotation_panel_20261007_v1';R=B/'multiscale_20261007_v1';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_multiscale_20261007_v1')
for p in [R/'scripts',R/'results',R/'figures',R/'logs',R/'standards',D/'graphs']:p.mkdir(parents=True,exist_ok=True)
src=B/'v3_pilot_20261006/vendor/quiche_graph.py';spec=importlib.util.spec_from_file_location('quiche_graph_native',src);qu=importlib.util.module_from_spec(spec);spec.loader.exec_module(qu)
oldcfg=json.loads((A/'standards/QUICHE_RETRY_PARAMETERS.json').read_text());assert hashlib.sha256(src.read_bytes()).hexdigest()==oldcfg['source_sha256']
ref=ad.read_h5ad('/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1/GSM9815454_panel_annotated_v1.h5ad',backed='r');allobs=ref.obs.copy();full=allobs.loc[allobs.qc_pass].copy();fullxy=ref.obsm['spatial'][allobs.qc_pass.to_numpy()];ref.file.close()
pilot=pd.read_csv(A/'results/pilot6000_annotations.tsv.gz',sep='\t',index_col=0);ids=pilot.index;ix=full.index.get_indexer(ids);assert (ix>=0).all();obs=full.loc[ids].copy();xy=fullxy[ix];labels=obs.annotation.astype(str).to_numpy();full_labels=full.annotation.astype(str).to_numpy();classes=sorted(set(full_labels));fullcodes=np.array([classes.index(x) for x in full_labels]);onehot=sparse.csr_matrix((np.ones(len(full)),(np.arange(len(full)),fullcodes)),shape=(len(full),len(classes)));hot=onehot[ix]
labelhash=hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(ids,labels)).encode()).hexdigest();assert labelhash==oldcfg['labels_sha256']
cfg=dict(sample='M242T/GSM9815454',n_biological_samples=1,focal_cells=6000,primary_neighbor_pool='same original 6000 cells',edge_sensitivity_neighbor_pool=f'all {len(full)} QC cells; same 6000 focal cells',labels_sha256=labelhash,source_sha256=oldcfg['source_sha256'],fixed_radius_um=[100,200,300],fixed_k=[15,30,60],exclude_self=True,min_neighbors=3,radius_definition='Euclidean distance < radius; no count cap',knn_definition='exact k other cells; no distance cap',radius_native='official compute_spatial_neighbors(n_neighbors=None, delaunay=False)',knn_native='official spatial_niches_khop(k=15/30/60,khop=1,radius=inf)',context_implementation='scipy cKDTree; independent brute-force Euclidean checks',uncertain_labels='retained as separate category',metrics='neighbor counts/extent, total variation, per-label rank correlation, focal lineage composition, window truncation sensitivity',inferential_tests=False,differential_abundance=False,CGNL1='not measured',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'standards/MULTISCALE_SPECIFICATION.json').write_text(json.dumps(cfg,indent=2))
container=ad.AnnData(X=sparse.csr_matrix((len(obs),1)),obs=obs[['sample','annotation']].copy());container.obs['annotation']=pd.Categorical(labels);container.obsm['spatial']=xy
P={};N={};stats=[];focal=[];validation=[];percell=pd.DataFrame(index=ids);percell['focal_label']=labels
tree=cKDTree(fullxy);auditix=np.sort(np.random.default_rng(20261007).choice(6000,20,replace=False))
def graph_from_rows(rows,ncols):
 lens=np.array([len(x) for x in rows]);return sparse.csr_matrix((np.ones(lens.sum(),dtype=np.float32),np.concatenate(rows),np.r_[0,np.cumsum(lens)]),shape=(6000,ncols))
def analyze(name,G,poolxy,poolhot,poollabels,selfindex,native_profile=None):
 assert G.shape==(6000,len(poolxy));G=G.tocsr();G.eliminate_zeros();assert np.all(np.asarray(G[np.arange(6000),selfindex]).ravel()==0)
 count=np.diff(G.indptr);assert np.all(count>=3);fractions=(G@poolhot).toarray()/count[:,None];assert np.isfinite(fractions).all() and np.allclose(fractions.sum(1),1)
 if native_profile is not None:assert np.allclose(fractions,native_profile.reindex(index=ids,columns=classes,fill_value=0).to_numpy())
 maximum=np.array([np.linalg.norm(poolxy[G.indices[G.indptr[i]:G.indptr[i+1]]]-xy[i],axis=1).max() for i in range(6000)])
 for i in auditix:
  d=np.linalg.norm(poolxy-xy[i],axis=1);d[selfindex[i]]=np.inf
  if 'radius' in name:j=np.flatnonzero(d<float(name.split('_')[-1]))
  else:j=np.argsort(d,kind='stable')[:int(name.split('_')[-1])]
  expected=pd.Series(poollabels[j]).value_counts(normalize=True).reindex(classes,fill_value=0).to_numpy();assert np.allclose(expected,fractions[i]),(name,i);assert len(j)==count[i]
 pd.DataFrame(fractions,index=ids,columns=classes).to_csv(R/f'results/{name}_composition.tsv.gz',sep='\t',index_label='cell_id');sparse.save_npz(D/'graphs'/f'{name}.npz',G)
 P[name]=fractions;N[name]=count;percell[name+'_neighbors']=count;percell[name+'_max_distance_um']=maximum
 row=dict(condition=name,cells=6000,neighbors_min=int(count.min()),neighbors_q25=float(np.quantile(count,.25)),neighbors_median=float(np.median(count)),neighbors_q75=float(np.quantile(count,.75)),neighbors_max=int(count.max()),distance_q25_um=float(np.quantile(maximum,.25)),distance_median_um=float(np.median(maximum)),distance_q75_um=float(np.quantile(maximum,.75)),distance_max_um=float(maximum.max()),mean_uncertain_fraction=float(fractions[:,classes.index('uncertain')].mean()),mean_hepatocyte_fraction=float(fractions[:,classes.index('hepatocyte_lineage')].mean()),endo_fraction_sd=float(fractions[:,classes.index('endothelial')].std()),mural_fraction_sd=float(fractions[:,classes.index('mural')].std()))
 stats.append(row);validation.append(dict(condition=name,finite_and_unit_row_sum=True,exclude_self=True,independent_20_cells_passed=True,native_profile_matches=native_profile is not None))
 baseline=np.bincount([classes.index(l) for l in poollabels],minlength=len(classes))/len(poollabels)
 for c in classes:
  mask=labels==c
  if mask.any():
   for j,n in enumerate(classes):focal.append(dict(condition=name,focal_label=c,focal_cells=int(mask.sum()),neighbor_label=n,mean_fraction=float(fractions[mask,j].mean()),pool_fraction=float(baseline[j]),difference_from_pool_percentage_points=float((fractions[mask,j].mean()-baseline[j])*100)))
 print('COMPLETE',json.dumps(row),flush=True)
for rad in [100,200,300]:
 a=container.copy();qu.compute_spatial_neighbors(a,radius=rad,n_neighbors=None,delaunay=False,fov_key='sample');G=a.obsp['spatial_connectivities'].tocsr();G.setdiag(0);G.eliminate_zeros();G.data[:]=1;analyze(f'window_radius_{rad}',G,xy,hot,labels,np.arange(6000))
for k in [15,30,60]:
 prof,nn=qu.spatial_niches_khop(container,radius=float('inf'),k=k,khop=1,labels_key='annotation',fov_key='sample',min_cell_threshold=3,n_jobs=8);G=graph_from_rows(nn['GSM9815454'],6000);analyze(f'window_knn_{k}',G,xy,hot,labels,np.arange(6000),prof)
baseline=pd.read_csv(A/'results/QUICHE_niches_radius100.tsv.gz',sep='\t',index_col=0).reindex(index=ids,columns=classes,fill_value=0);assert np.allclose(P['window_knn_30'],baseline.to_numpy());print('ORIGINAL_K30_EXACTLY_REPRODUCED',flush=True)
for rad in [100,200,300]:
 rows=tree.query_ball_point(xy,r=rad,workers=8);rows=[np.asarray(row,dtype=np.int32) for row in rows]
 rows=[row[(row!=ix[i])&(np.linalg.norm(fullxy[row]-xy[i],axis=1)<rad)] for i,row in enumerate(rows)];G=graph_from_rows(rows,len(full));analyze(f'context_radius_{rad}',G,fullxy,onehot,full_labels,ix)
distance,neighbors=tree.query(xy,k=61,workers=8);rows=[nn[nn!=ix[i]][:60].astype(np.int32) for i,nn in enumerate(neighbors)]
for k in [15,30,60]:
 G=graph_from_rows([row[:k] for row in rows],len(full));analyze(f'context_knn_{k}',G,fullxy,onehot,full_labels,ix)
pairwise=[];correlations=[];edges=[]
for pool in ['window','context']:
 for family,values in [('radius',[100,200,300]),('knn',[15,30,60])]:
  for v1,v2 in itertools.combinations(values,2):
   a=f'{pool}_{family}_{v1}';b=f'{pool}_{family}_{v2}';tv=.5*np.abs(P[a]-P[b]).sum(1);pairwise.append(dict(pool=pool,family=family,scale1=v1,scale2=v2,mean_TV=float(tv.mean()),median_TV=float(np.median(tv)),p90_TV=float(np.quantile(tv,.9)),fraction_TV_above_0_1=float(np.mean(tv>.1)),exact_equal_fraction=float(np.mean(tv<1e-12))));percell[f'{pool}_{family}_{v1}_vs_{v2}_TV']=tv
   for j,c in enumerate(classes):
    rho=float(spearmanr(P[a][:,j],P[b][:,j]).statistic) if np.ptp(P[a][:,j]) and np.ptp(P[b][:,j]) else None
    correlations.append(dict(pool=pool,family=family,scale1=v1,scale2=v2,label=c,spearman_rho=rho,mean_absolute_change_percentage_points=float(np.abs(P[a][:,j]-P[b][:,j]).mean()*100)))
for family,values in [('radius',[100,200,300]),('knn',[15,30,60])]:
 for v in values:
  a=f'window_{family}_{v}';b=f'context_{family}_{v}';tv=.5*np.abs(P[a]-P[b]).sum(1);edges.append(dict(family=family,scale=v,mean_TV=float(tv.mean()),median_TV=float(np.median(tv)),p90_TV=float(np.quantile(tv,.9)),fraction_composition_changed=float(np.mean(tv>1e-12)),mean_extra_neighbors=float((N[b]-N[a]).mean()),median_extra_neighbors=float(np.median(N[b]-N[a])),fraction_neighbor_count_changed=float(np.mean(N[b]!=N[a]))));percell[f'{family}_{v}_window_vs_context_TV']=tv
outside_idx=np.flatnonzero(~np.isin(np.arange(len(full)),ix));edge_distance=cKDTree(fullxy[outside_idx]).query(xy,k=1,workers=8)[0];percell['nearest_excluded_window_cell_um']=edge_distance;interior=edge_distance>=300;assert interior.any();core=[]
for row in pairwise:
 a=f"{row['pool']}_{row['family']}_{row['scale1']}";b=f"{row['pool']}_{row['family']}_{row['scale2']}";tv=.5*np.abs(P[a]-P[b]).sum(1);core.append(dict(**row,interior_cells=int(interior.sum()),interior_mean_TV=float(tv[interior].mean()),interior_median_TV=float(np.median(tv[interior]))))
for rad in [100,200,300]:assert np.allclose(P[f'window_radius_{rad}'][interior],P[f'context_radius_{rad}'][interior])
for name,data in [('scale_summary',stats),('scale_comparisons',pairwise),('composition_correlations',correlations),('window_edge_sensitivity',edges),('interior_sensitivity',core),('focal_lineage_composition',focal)]:pd.DataFrame(data).to_csv(R/f'results/{name}.tsv',sep='\t',index=False)
percell.to_csv(R/'results/per_cell_multiscale_metrics.tsv.gz',sep='\t',index_label='cell_id');(R/'results/NUMERICAL_VALIDATION.json').write_text(json.dumps(validation,indent=2))
out=dict(exit_code=0,focal_cells=6000,context_cells=len(full),conditions=len(stats),interior_300um_cells=int(interior.sum()),original_k30_exactly_reproduced=True,all_conditions_independently_validated=True,labels_sha256=labelhash,no_inferential_tests=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(R/'results/MULTISCALE_COMPLETE.json').write_text(json.dumps(out,indent=2));print('ALL_COMPLETE',json.dumps(out),flush=True)
