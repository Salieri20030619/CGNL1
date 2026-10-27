import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['NUMBA_NUM_THREADS']='2'
import pathlib,json,numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2/pair_validation');C=O/'epithelial_states';C.mkdir(exist_ok=True);p=json.loads((O/'PLAN.json').read_text());heldout=sorted(set(p['genes']+sum(p['programs'].values(),[])+p['gastric_proxy']));plan={'trigger':'Independent broad epithelial OPLAH direction contradicts bulk; distinguish subtype mixture from within-state expression.','heldout':heldout,'HVG':2000,'batch_key':'sample','PCA':30,'neighbors':30,'leiden':0.5,'seed':23,'integration':'No batch correction that might erase PT/PM; report sample dominance.','shared_state_rule':'At least2 PT and2 PM samples with>=20 cells each; sensitivity>=50. Descriptive only.','malignancy':'No CNV confirmation; do not call epithelial clusters malignant.'};(C/'PLAN.json').write_text(json.dumps(plan,indent=2))
parts=[];metas=[];ref=None
for f in sorted(O.glob('*_epithelial_counts.npz')):
 pre=f.name[:-len('_epithelial_counts.npz')];ft=pd.read_csv(O/(pre+'_features.tsv.gz'),sep='\t');m=pd.read_csv(O/(pre+'_epithelial_metadata.tsv.gz'),sep='\t');x=sparse.load_npz(f).T.tocsr();assert len(m)==x.shape[0]
 if ref is None:ref=ft;assert ref.gene_id.is_unique
 else:assert np.array_equal(ref.gene_id.values,ft.gene_id.values)
 m.index=pre+'__'+m.barcode.astype(str);parts.append(x);metas.append(m)
counts=sparse.vstack(parts,format='csr').astype(np.float32);obs=pd.concat(metas);A=ad.AnnData(counts.copy(),obs=obs,var=ref.set_index('gene_id'));sc.settings.n_jobs=2;sc.pp.normalize_total(A,target_sum=1e4);sc.pp.log1p(A);sc.pp.highly_variable_genes(A,n_top_genes=2000,flavor='seurat',batch_key='sample');keep=A.var.highly_variable.values&~A.var.symbol.isin(heldout).values&~A.var.symbol.str.startswith(('MT-','RPS','RPL')).values;A.var.to_csv(C/'HVG_audit.tsv',sep='\t');assert not set(A.var.loc[keep,'symbol'])&set(heldout)
Z=A[:,keep].copy();sc.pp.scale(Z,max_value=10);sc.tl.pca(Z,n_comps=30,svd_solver='arpack',random_state=23);sc.pp.neighbors(Z,n_neighbors=30,n_pcs=30,random_state=23);sc.tl.leiden(Z,resolution=.5,random_state=23,key_added='cluster');obs['cluster']=Z.obs.cluster.values;obs=obs.reset_index(names='cell_id');obs.to_csv(C/'cell_clusters.tsv.gz',sep='\t',index=False)
ct=pd.crosstab(obs['sample'],obs.cluster);ct.to_csv(C/'sample_cluster_counts.tsv',sep='\t');qq=[]
for k in ct:
 v=ct[k];qq.append({'cluster':str(k),'cells':int(v.sum()),'max_sample_fraction':float(v.max()/v.sum()),'PT_samples_ge20':int(((v.index.str.contains('_Ca_'))&(v>=20)).sum()),'PM_samples_ge20':int(((v.index.str.contains('_F_'))&(v>=20)).sum())})
qc=pd.DataFrame(qq);qc.to_csv(C/'cluster_QC.tsv',sep='\t',index=False)
markers=['EPCAM','KRT8','KRT18','KRT19','KRT7','MSLN','CALB2','WT1','UPK3B','PTPRC','LST1','COL1A1','VIM','CDH1','MUC1','MUC5AC','MUC6','PGC','LIPF','GKN1','GKN2','CLDN18','CEACAM5','CEACAM6','TFF1','TFF2','TFF3','MKI67','TOP2A'];wanted=sorted(set(markers+heldout));sy=ref.symbol.values;rows=[]
for (sample,k),ix in obs.groupby(['sample','cluster'],observed=True).groups.items():
 ix=np.array(list(ix));pb=np.asarray(counts[ix].sum(axis=0)).ravel();lib=pb.sum();n=len(ix)
 for g in wanted:
  gi=np.flatnonzero(sy==g)
  if not len(gi):continue
  cellcounts=np.asarray(counts[ix][:,gi].sum(axis=1)).ravel();r={'sample':sample,'group':obs.loc[ix[0],'group'],'cluster':str(k),'gene':g,'n_cells':n,'log2CPM':float(np.log2(1+pb[gi].sum()/max(lib,1)*1e6)),'detection_fraction':float((cellcounts>0).mean()),'mean_CP10K':float(np.mean(cellcounts/obs.loc[ix,'umi'].values*1e4))};rows.append(r)
a=pd.DataFrame(rows);a.to_csv(C/'state_sample_expression.tsv.gz',sep='\t',index=False);summary=[]
for cutoff in [20,50]:
 for k in sorted(a.cluster.unique()):
  for g in p['genes']:
   z=a[(a.cluster==k)&(a.gene==g)&(a.n_cells>=cutoff)];pt=z[z.group=='PT'];pm=z[z.group=='PM'];d=pm.log2CPM.values[:,None]-pt.log2CPM.values[None,:];summary.append({'cluster':k,'gene':g,'min_cells':cutoff,'PT_samples':len(pt),'PM_samples':len(pm),'eligible_shared':len(pt)>=2 and len(pm)>=2,'PT_median':pt.log2CPM.median(),'PM_median':pm.log2CPM.median(),'delta_log2CPM':pm.log2CPM.median()-pt.log2CPM.median(),'delta_mean_CP10K':pm.mean_CP10K.median()-pt.mean_CP10K.median(),'cross_comparisons_lower':int((d<0).sum()),'cross_comparisons_total':d.size})
s=pd.DataFrame(summary);s.to_csv(C/'within_state_candidate_comparison.tsv',sep='\t',index=False);(C/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'cells':len(obs),'samples':obs['sample'].nunique(),'clusters':obs.cluster.nunique(),'HVG_after_holdout':int(keep.sum()),'eligible_shared_states_min20':int(s[(s.min_cells==20)&s.eligible_shared].cluster.nunique()),'scope':'descriptive sample-level, not malignancy-confirmed, no P/FDR'},indent=2));print(qc.to_string(index=False));print(s[s.eligible_shared].to_string(index=False))
