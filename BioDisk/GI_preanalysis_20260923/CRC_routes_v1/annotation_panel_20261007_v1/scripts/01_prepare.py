from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
import json,gzip,io,hashlib,datetime
import numpy as np,pandas as pd,h5py,anndata as ad,scanpy as sc
from scipy import sparse
from sklearn.metrics import balanced_accuracy_score
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'annotation_panel_20261007_v1';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1')
oldD=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
spec=json.loads((B/'identity_review_target_blind_20261002_v1/SPECIFICATION.json').read_text())
excluded=set(spec['excluded_target_panel'])|set(spec['extra_exclusion_sensitivity'])
panels={
 'hepatocyte':('HNF1A G6PC CYP2E1 CYP2C9 FABP1 APOB FGB LECT2','APOC2 CFHR5 PLG ORM1 GLUL','EPCAM KRT19 PECAM1 CD3D CD79A'),
 'endothelial':('PECAM1 VWF KDR FLT1 CLEC14A PTPRB ROBO4','PLVAP CD34 EPAS1 ACKR1 CLEC4G LYVE1 FCN3','CD3D CD79A EPCAM'),
 'mural':('RGS5 MYH11 ACTA2 TAGLN ITGA8','NOTCH3 PLN THSD4','PECAM1 VWF EPCAM CD3D'),
 'fibroblast':('COL1A1 COL1A2 DCN PDGFRA CLEC3B','FAP POSTN BGN MMP2 THBS2 HGF FOXF1','PECAM1 VWF EPCAM CD3D'),
 'epithelial':('EPCAM KRT19 KRT7 TACSTD2 CFTR','SOX9 REG3A MUC6 LGR5 PROM1 EGFR ERBB2','PECAM1 CD3D CD79A'),
 'T':('CD3D CD3E CD7 THEMIS CD8B','CD4 CD8A IL7R TCF7 TRDC IL32','MS4A1 CD79A TPSB2'),
 'NK':('NKG7 GNLY PRF1 FGFBP2 KLRC1','EOMES GZMB GZMH XCL1 XCL2 NCAM1','CD3D CD3E CD79A'),
 'B':('MS4A1 CD19 CD79A CD79B BANK1','FCRL1 CR2 TNFRSF13C IGHM IGHD','CD3D MZB1 TPSB2'),
 'plasma':('MZB1 TNFRSF17 SDC1','CD38 IGHG3 IGHM','CD3D MS4A1'),
 'myeloid':('CD68 CD163 C1QA MERTK FCGR1A CD14','LYZ MRC1 MARCO FOLR2 CD5L TIMD4 STAB1 TREM2 FCGR3A FCN1 ITGAM','CD3D CD79A EPCAM'),
 'dendritic':('CD1C FCER1A CLEC9A CLEC4C LILRA4','CD1E FLT3 LAMP3 ITGAX CD209 HLA-DRA','CD3D CD79A'),
 'mast':('TPSB2 KIT','','CD3D CD79A'),
 'neutrophil':('FCGR3B CSF3R CEACAM8 CXCR2','S100A9 LTF FUT4','CD3D CD79A')}
with gzip.open(oldD/'downloads/GSM9815454_cells.parquet.gz','rb') as h:obs=pd.read_parquet(io.BytesIO(h.read())).set_index('cell_id')
raw=B/'spatial_v3/raw/GSM9815454_tumor1_M242T_cell_feature_matrix.h5'
with h5py.File(raw) as f:
 g=f['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr();names=np.array([v.decode() for v in g['features/name'][:]]);types=np.array([v.decode() for v in g['features/feature_type'][:]]);ids=[v.decode() for v in g['barcodes'][:]]
obs=obs.loc[ids].copy();keep=types=='Gene Expression';X=X[:,keep];names=names[keep]
assert np.array_equal(np.asarray(X.sum(1)).ravel(),obs.transcript_counts.to_numpy())
assert obs.index.is_unique and len(names)==462
obs['qc_pass']=(obs.transcript_counts>=50)&(obs.nucleus_count==1)
obs['negative_control_fraction']=(obs.control_probe_counts+obs.control_codeword_counts)/np.maximum(obs.total_counts,1)
obs['area_flag']=(obs.nucleus_area>obs.cell_area)|(obs.cell_area<=0)
obs['qc_reason']=np.where(obs.qc_pass,'pass',np.where(obs.nucleus_count!=1,'nucleus_count_not_1','transcripts_below_50'))
old=ad.read_h5ad(oldD/'inputs/GSM9815454_quiche_technical.h5ad',backed='r')
obs['old_pilot_label']='not_in_pilot';obs.loc[old.obs_names,'old_pilot_label']=old.obs.candidate_identity.astype(str).to_numpy();old.file.close()
a=ad.AnnData(X.astype(np.float32),obs=obs,var=pd.DataFrame(index=names));a.obsm['spatial']=obs[['x_centroid','y_centroid']].to_numpy();a.write_h5ad(D/'inputs/full_counts.h5ad',compression='gzip')
rows=[]
for label,groups in panels.items():
 for role,genes in zip(['core','support','conflict'],groups):
  for gene in genes.split():rows.append(dict(label=label,role=role,gene=gene,measured=gene in names,allowed_identity=gene in names and gene not in excluded))
pd.DataFrame(rows).to_csv(R/'standards/panel_marker_standard.tsv',sep='\t',index=False)
cfg=dict(sample='GSM9815454_M242T',biological_samples=1,panels=panels,excluded_targets=sorted(excluded),normalization='log1p 10000 / Gene Expression counts only',qc='transcript_counts>=50 and nucleus_count==1; all original cells retained in export',clustering='expression-only scaled PCA30 cosine k30 Leiden 0.3/0.6/1.0, seed20261007',annotation='review expression clusters, core-marker enrichment and independent reference agreement; conflicts remain uncertain; no malignancy assignment',reference='donor-balanced Pearson centroids on measured identity core/support union; no re-normalization of already log-normalized references',rare_lineages='review only if at least two independent core markers enriched; T and NK may remain T_NK; DC remains myeloid if unsupported',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'standards/ANNOTATION_SPECIFICATION.json').write_text(json.dumps(cfg,indent=2))
q=a[a.obs.qc_pass].copy();sc.pp.normalize_total(q,target_sum=1e4);sc.pp.log1p(q)
eligible=~q.var_names.isin(excluded);z=q[:,eligible].copy();sc.pp.filter_genes(z,min_cells=10);sc.pp.scale(z,max_value=10);sc.tl.pca(z,n_comps=30,random_state=20261007,svd_solver='arpack');sc.pp.neighbors(z,n_neighbors=30,n_pcs=30,metric='cosine',random_state=20261007)
for res,key in [(.3,'cluster03'),(.6,'cluster06'),(1.,'cluster10')]:sc.tl.leiden(z,resolution=res,key_added=key,random_state=20261007,n_iterations=2,flavor='igraph');q.obs[key]=z.obs[key].to_numpy()
q.obsm['X_pca']=z.obsm['X_pca'];q.obsp['connectivities']=z.obsp['connectivities'];q.obsp['distances']=z.obsp['distances']
pd.Series(z.var_names,name='gene').to_csv(R/'standards/clustering_genes.tsv',sep='\t',index=False)
del z
print('CLUSTERS',q.obs.cluster06.value_counts().to_dict(),flush=True)
def corr(Q,P):
 Q=np.asarray(Q,dtype=np.float64);P=np.asarray(P,dtype=np.float64);Q=Q-Q.mean(1,keepdims=True);P=P-P.mean(1,keepdims=True);return Q@P.T/np.maximum(np.linalg.norm(Q,axis=1)[:,None]*np.linalg.norm(P,axis=1)[None,:],1e-12)
marker_genes=sorted({g for groups in panels.values() for gs in groups[:2] for g in gs.split()}&set(names)-excluded)
for refname in ['normal_liver','multiome']:
 ref=ad.read_h5ad(B/f'identity_review_target_blind_20261002_v1/references/{refname}_target_excluded.h5ad',backed='r');gs=[g for g in marker_genes if g in ref.var_names];m=ref.obs.copy();V=ref[:,gs].X;V=V.toarray() if sparse.issparse(V) else np.asarray(V);ref.file.close()
 profiles=[];pkeys=[]
 for (donor,label),ix in m.groupby(['donor','label'],observed=True).indices.items():profiles.append(V[ix].mean(0));pkeys.append((str(donor),str(label)))
 profiles=np.array(profiles);pk=pd.DataFrame(pkeys,columns=['donor','label'])
 def centroids(ex=None):
  good=pk.donor.ne(ex);labels=sorted(pk.loc[good,'label'].unique());return labels,np.array([profiles[good&pk.label.eq(l)].mean(0) for l in labels])
 val=[];predrows=[]
 for donor in sorted(m.donor.unique()):
  labs,P=centroids(donor);mask=m.donor.eq(donor).to_numpy();pred=np.array(labs)[corr(V[mask],P).argmax(1)];truth=m.loc[mask,'label'].astype(str).to_numpy();val.append(dict(reference=refname,donor=donor,cells=int(mask.sum()),accuracy=float(np.mean(pred==truth)),balanced_accuracy=float(balanced_accuracy_score(truth,pred))));predrows.extend(zip([donor]*len(pred),truth,pred))
 pd.DataFrame(val).to_csv(R/f'results/{refname}_panel_LODO.tsv',sep='\t',index=False)
 pv=pd.DataFrame(predrows,columns=['donor','truth','prediction']);pd.crosstab(pv.truth,pv.prediction).to_csv(R/f'results/{refname}_panel_LODO_confusion.tsv',sep='\t')
 labs,P=centroids();pd.DataFrame(P,index=labs,columns=gs).to_csv(R/f'standards/{refname}_panel_centroids.tsv',sep='\t');Q=q[:,gs].X.toarray();S=corr(Q,P);order=np.argsort(S,axis=1);q.obs[refname+'_prediction']=np.array(labs)[order[:,-1]];q.obs[refname+'_score']=S[np.arange(len(S)),order[:,-1]];q.obs[refname+'_margin']=S[np.arange(len(S)),order[:,-1]]-S[np.arange(len(S)),order[:,-2]]
 cr=[]
 for cl,ix in q.obs.groupby('cluster06',observed=True).indices.items():
  ss=corr(Q[ix].mean(0)[None,:],P)[0];o=np.argsort(ss);cr.append(dict(cluster=str(cl),prediction=labs[o[-1]],score=ss[o[-1]],margin=ss[o[-1]]-ss[o[-2]],second=labs[o[-2]]))
 pd.DataFrame(cr).to_csv(R/f'results/{refname}_cluster_reference.tsv',sep='\t',index=False)
 pd.crosstab(q.obs.cluster06,q.obs[refname+'_prediction'],normalize='index').to_csv(R/f'results/{refname}_cluster_reference_fractions.tsv',sep='\t')
 print('REFERENCE_DONE',refname,len(gs),flush=True)
means=[];fracs=[];markers=[];globalmean=np.asarray(q.X.mean(0)).ravel();ng=q.n_obs
for cl,ix in q.obs.groupby('cluster06',observed=True).indices.items():
 Z=q.X[ix];mean=np.asarray(Z.mean(0)).ravel();pct=Z.getnnz(axis=0)/len(ix);other=(globalmean*ng-mean*len(ix))/(ng-len(ix));delta=mean-other;means.append(pd.Series(mean,index=names,name=str(cl)));fracs.append(pd.Series(pct,index=names,name=str(cl)))
 for j in np.argsort(delta)[::-1]:markers.append(dict(cluster=str(cl),cells=len(ix),gene=names[j],mean_log=mean[j],delta_log_mean=delta[j],fraction_detected=pct[j],identity_allowed=names[j] not in excluded))
pd.DataFrame(means).to_csv(R/'results/cluster_mean_log.tsv',sep='\t');pd.DataFrame(fracs).to_csv(R/'results/cluster_detection_fraction.tsv',sep='\t');pd.DataFrame(markers).to_csv(R/'results/cluster_markers_descriptive.tsv',sep='\t',index=False)
q.write_h5ad(D/'cache/clustered_log_normalized.h5ad',compression='gzip')
q.obs.to_csv(R/'results/preannotation_cell_evidence.tsv.gz',sep='\t',index_label='cell_id')
out=dict(exit_code=0,all_cells=a.n_obs,qc_cells=q.n_obs,genes=len(names),clusters=q.obs.cluster06.nunique(),raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),excluded=sorted(excluded),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'results/PREPARE_COMPLETE.json').write_text(json.dumps(out,indent=2));print('PREPARE_COMPLETE',json.dumps(out),flush=True)
