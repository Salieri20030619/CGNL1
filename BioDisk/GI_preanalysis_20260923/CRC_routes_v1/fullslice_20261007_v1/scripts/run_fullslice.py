from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
import json,hashlib,datetime,importlib.metadata as md
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');A=B/'annotation_panel_20261007_v1';M=B/'multiscale_20261007_v1';R=B/'fullslice_20261007_v1';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_fullslice_20261007_v1');D.mkdir(exist_ok=True);(D/'graphs').mkdir(exist_ok=True)
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
spec=dict(version='20261007_v1',scope='full M242T slice, previously reviewed labels; descriptive neighborhood analysis',n_independent_patients=1,primary='knn_30',sensitivity=['knn_15','knn_60'],physical_distance_complement='radius_100',neighbor_pool='all QC cells in same slice',exclude_self=True,min_neighbors=3,labels='frozen panel-adapted reviewed annotations; uncertain retained',input='/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1/GSM9815454_panel_annotated_v1.h5ad',engine='vectorized cKDTree neighborhood composition; definition validated against QUICHE native module in prior run; not full QUICHE differential abundance',grouping=dict(method='KMeans on square-root cell-type fractions',n_groups=6,n_init=10,primary_seed=20261007,replicate_seed=20261008,sensitivity='project changed-scale features to frozen k30 centroids',interpretation='six operational descriptive groups; no optimal-K claim or pathology assignment'),statistics='descriptive cell-level effects; no biological P or FDR with n=1',CGNL1='not measured; no direct target validation',prior_pilot_overlap='same sample, not independent replication',download_policy='direct_only; no new external data required',frozen_utc=now())
(R/'standards/FULLSLICE_SPECIFICATION.json').write_text(json.dumps(spec,indent=2))
ref=ad.read_h5ad(spec['input'],backed='r');obs=ref.obs.loc[ref.obs.qc_pass].copy();xy=ref.obsm['spatial'][ref.obs.qc_pass.to_numpy()];ref.file.close();n=len(obs);assert n==75094 and obs.index.is_unique and np.isfinite(xy).all()
labels=obs.annotation.astype(str).to_numpy();classes=sorted(set(labels));codes=pd.Categorical(labels,categories=classes).codes;hot=sparse.csr_matrix((np.ones(n),(np.arange(n),codes)),shape=(n,len(classes)));tree=cKDTree(xy)
labelhash=hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(obs.index,labels)).encode()).hexdigest();(R/'standards/FULL_LABEL_FREEZE.json').write_text(json.dumps(dict(cells=n,labels_sha256=labelhash,utc=now()),indent=2))
pilot=pd.read_csv(A/'results/pilot6000_annotations.tsv.gz',sep='\t',index_col=0);pi=obs.index.get_indexer(pilot.index);assert (pi>=0).all();oldcfg=json.loads((M/'standards/MULTISCALE_SPECIFICATION.json').read_text());assert hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(pilot.index,labels[pi])).encode()).hexdigest()==oldcfg['labels_sha256']
audits=np.random.default_rng(20261007).choice(n,30,replace=False);summaries=[];validations=[];focal=[];P={};meta=obs[['annotation','cluster_identity','review_cluster','annotation_confidence','mixed_marker_flag','transcript_counts','negative_control_fraction','segmentation_method','cell_area','nucleus_area','x_centroid','y_centroid']].copy()
def fromrows(rows):
 lengths=np.array([len(v) for v in rows]);return sparse.csr_matrix((np.ones(lengths.sum(),dtype=np.float32),np.concatenate(rows),np.r_[0,np.cumsum(lengths)]),shape=(n,n))
def process(name,G):
 G=G.tocsr();G.eliminate_zeros();assert np.all(G.diagonal()==0);cnt=np.diff(G.indptr);keep=cnt>=3;frac=(G@hot).toarray()/np.maximum(cnt,1)[:,None];frac[~keep]=np.nan;assert np.isfinite(frac[keep]).all() and np.allclose(frac[keep].sum(1),1)
 # Direct independent checks use all point distances, never reuse the search result.
 for i in audits:
  dd=np.linalg.norm(xy-xy[i],axis=1);dd[i]=np.inf
  if name.startswith('knn'):jj=np.argsort(dd,kind='stable')[:int(name.split('_')[-1])]
  else:jj=np.flatnonzero(dd<100)
  assert len(jj)==cnt[i]
  if keep[i]:assert np.allclose(frac[i],np.bincount(codes[jj],minlength=len(classes))/len(jj))
 old=pd.read_csv(M/f'results/context_{name}_composition.tsv.gz',sep='\t',index_col=0).reindex(index=pilot.index,columns=classes,fill_value=0);assert np.allclose(frac[pi],old.to_numpy())
 maximum=np.zeros(n)
 for start in range(0,n,2048):
  for i in range(start,min(start+2048,n)):
   jj=G.indices[G.indptr[i]:G.indptr[i+1]];maximum[i]=np.linalg.norm(xy[jj]-xy[i],axis=1).max() if len(jj) else np.nan
 sparse.save_npz(D/'graphs'/f'{name}.npz',G);pd.DataFrame(frac,index=obs.index,columns=classes).to_csv(R/f'results/{name}_composition.tsv.gz',sep='\t',index_label='cell_id');P[name]=frac;meta[name+'_neighbors']=cnt;meta[name+'_max_distance_um']=maximum;meta[name+'_valid']=keep
 row=dict(condition=name,centers=n,retained_centers=int(keep.sum()),neighbors_min=int(cnt.min()),neighbors_q25=float(np.quantile(cnt,.25)),neighbors_median=float(np.median(cnt)),neighbors_q75=float(np.quantile(cnt,.75)),neighbors_max=int(cnt.max()),distance_median_um=float(np.nanmedian(maximum)),distance_p95_um=float(np.nanquantile(maximum,.95)),mean_uncertain_fraction=float(np.nanmean(frac[:,classes.index('uncertain')])));summaries.append(row)
 overall=np.nanmean(frac,axis=0)
 for c in classes:
  mask=(labels==c)&keep
  if not mask.any():continue
  for j,target in enumerate(classes):focal.append(dict(condition=name,focal_label=c,focal_cells=int(mask.sum()),neighbor_label=target,mean_fraction=float(frac[mask,j].mean()),all_center_neighborhood_mean=float(overall[j]),difference_percentage_points=float(100*(frac[mask,j].mean()-overall[j]))))
 validations.append(dict(condition=name,independent_30_rows_passed=True,all_6000_previous_context_rows_reproduced=True,self_excluded=True,valid_rows_finite_unit_sum=True));print('NEIGHBORHOOD_COMPLETE',json.dumps(row),flush=True)
dist,knn=tree.query(xy,k=61,workers=8);rows=[nn[nn!=i][:60].astype(np.int32) for i,nn in enumerate(knn)]
for k in [30,15,60]:process(f'knn_{k}',fromrows([v[:k] for v in rows]))
rows=tree.query_ball_point(xy,100,workers=8);rows=[np.asarray(v,dtype=np.int32) for v in rows];rows=[v[(v!=i)&(np.linalg.norm(xy[v]-xy[i],axis=1)<100)] for i,v in enumerate(rows)];process('radius_100',fromrows(rows));del rows
z=np.sqrt(P['knn_30']);assert np.isfinite(z).all();km=KMeans(n_clusters=6,n_init=10,random_state=20261007);groups=km.fit_predict(z);alt=KMeans(n_clusters=6,n_init=10,random_state=20261008).fit_predict(z);seedari=float(adjusted_rand_score(groups,alt));cm=pd.crosstab(groups,alt).to_numpy();rr,cc=linear_sum_assignment(-cm);seedagree=float(cm[rr,cc].sum()/n)
meta['primary_niche']=['N'+str(v+1) for v in groups];meta['alternate_seed_cluster']=alt.astype(str);pd.DataFrame(km.cluster_centers_,index=['N'+str(i+1) for i in range(6)],columns=classes).to_csv(R/'standards/frozen_k30_sqrt_centroids.tsv',sep='\t');groupmeans=[];qc=[]
for i in range(6):
 mask=groups==i;mean=P['knn_30'][mask].mean(0);globalmean=P['knn_30'].mean(0);extra=mean-globalmean;best=np.argsort(extra)[::-1][:2];descriptor=' / '.join(classes[j] for j in best if extra[j]>.02) or 'near_average';groupmeans.append(dict(niche='N'+str(i+1),centers=int(mask.sum()),fraction=float(mask.mean()),descriptor=descriptor,**{c:float(v) for c,v in zip(classes,mean)}));qc.append(dict(niche='N'+str(i+1),transcripts_median=float(obs.iloc[np.flatnonzero(mask)].transcript_counts.median()),focal_uncertain_fraction=float(np.mean(labels[mask]=='uncertain')),negative_control_fraction_median=float(obs.iloc[np.flatnonzero(mask)].negative_control_fraction.median())))
stability=[]
for name in ['knn_15','knn_60','radius_100']:
 good=np.isfinite(P[name]).all(1);projection=np.full(n,-1,int);projection[good]=km.predict(np.sqrt(P[name][good]));meta[name+'_projected_niche']=np.where(good,['N'+str(v+1) for v in projection],'insufficient_neighbors');tv=.5*np.abs(P[name]-P['knn_30']).sum(1);meta[name+'_vs_k30_TV']=tv;stability.append(dict(condition=name,common_centers=int(good.sum()),same_frozen_group_fraction=float(np.mean(groups[good]==projection[good])),ARI_to_k30=float(adjusted_rand_score(groups[good],projection[good])),mean_TV=float(np.nanmean(tv)),median_TV=float(np.nanmedian(tv)),p90_TV=float(np.nanquantile(tv,.9))));pd.crosstab(meta.primary_niche,meta[name+'_projected_niche']).to_csv(R/f'results/{name}_niche_transitions.tsv',sep='\t')
meta.to_csv(R/'results/fullslice_cell_neighborhoods.tsv.gz',sep='\t',index_label='cell_id');pd.DataFrame(groupmeans).to_csv(R/'results/primary_niche_summary.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(R/'results/niche_QC_summary.tsv',sep='\t',index=False);pd.crosstab(meta.primary_niche,meta.annotation).to_csv(R/'results/niche_focal_identity_counts.tsv',sep='\t');pd.crosstab(meta.primary_niche,meta.segmentation_method).to_csv(R/'results/niche_segmentation_counts.tsv',sep='\t')
for name,data in [('scale_summary',summaries),('scale_stability',stability),('focal_lineage_composition',focal)]:pd.DataFrame(data).to_csv(R/f'results/{name}.tsv',sep='\t',index=False)
(R/'results/NUMERICAL_VALIDATION.json').write_text(json.dumps(validations,indent=2));(R/'standards/ENVIRONMENT.json').write_text(json.dumps({k:md.version(k) for k in ['numpy','pandas','scipy','scikit-learn','anndata']},indent=2))
embedding=ad.AnnData(X=P['knn_30'].astype(np.float32),obs=meta.copy(),var=pd.DataFrame(index=classes));embedding.obsm['spatial']=xy;embedding.uns['specification_json']=json.dumps(spec);embedding.write_h5ad(D/'M242T_fullslice_neighborhoods.h5ad',compression='gzip')
out=dict(exit_code=0,status='full_slice_descriptive_complete',sample='M242T',centers=n,conditions=4,groups=6,seed_ARI=seedari,seed_hungarian_agreement=seedagree,seed_stability_not_accuracy=True,all_previous_pilot_profiles_reproduced=True,full_label_sha256=labelhash,n_biological_samples=1,differential_abundance=False,CGNL1_direct_test=False,completed_utc=now());(R/'results/FULLSLICE_COMPLETE.json').write_text(json.dumps(out,indent=2));print('FULLSLICE_COMPLETE',json.dumps(out),flush=True)
