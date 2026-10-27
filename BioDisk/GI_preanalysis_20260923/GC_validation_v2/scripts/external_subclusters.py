import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['NUMBA_NUM_THREADS']='2'
import pathlib,json,tarfile,gzip,numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import io,sparse,stats
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_routes';C=O/'external_subclusters';C.mkdir(exist_ok=True);plan=json.loads((O/'PLAN.json').read_text());heldout=plan['genes']+plan['programs']['phagolysosomal'];(C/'PLAN.json').write_text(json.dumps({'trigger':'Distinguish within-lineage associations from macrophage subtype mixture','population':'existing C1Q macrophage-enriched gates; no candidate-based selection','heldout_from_clustering':heldout,'HVG':2000,'HVG_batch_key':'patient','normalization':'counts per10000 then log1p','PCA':30,'neighbors':30,'leiden_resolution':0.5,'seed':23,'integration':'no forced batch correction; donor dominance audited','association':'same two phagolysosomal links within each cluster, min20 cells per patient and min10 patients; BH across2*K including unavailable tests as1','composition_sensitivity':'all18 gene-program correlations in pooled macrophages, prior covariates plus<=5 PCs explaining>=80% of CLR cluster composition variance; BH within18'},indent=2))
parts=[];metas=[];refids=None
with tarfile.open(V/'02_raw/GSE228598_RAW.tar') as tar:
 for prefix in sorted(n[:-len('_matrix.mtx.gz')] for n in tar.getnames() if n.endswith('_matrix.mtx.gz')):
  def f(n):return gzip.GzipFile(fileobj=tar.extractfile(n))
  ft=pd.read_csv(f(prefix+'_features.tsv.gz'),sep='\t',header=None);ids=ft.iloc[:,0].astype(str).values;sy=ft.iloc[:,1].astype(str).values;bc=pd.read_csv(f(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0]
  if refids is None:refids=ids;symbols=sy
  else:assert np.array_equal(refids,ids)
  with threadpool_limits(limits=2):X=io.mmread(f(prefix+'_matrix.mtx.gz')).tocsr()
  m=pd.read_csv(V/'03_processed'/(prefix+'_gating_refined.tsv.gz'),sep='\t');m=m[m.label=='Macrophage_enriched'].copy();ix=pd.Index(bc).get_indexer(m.barcode);assert (ix>=0).all();parts.append(X[:,ix].T.tocsr());m.index=prefix+'__'+m.barcode.astype(str);metas.append(m)
  print('read',prefix,len(m),flush=True)
counts=sparse.vstack(parts,format='csr').astype(np.float32);obs=pd.concat(metas);assert counts.shape[0]==43113;v=pd.DataFrame({'symbol':symbols},index=refids);adata=ad.AnnData(counts.copy(),obs=obs,var=v);sc.settings.n_jobs=2;sc.pp.normalize_total(adata,target_sum=1e4);sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,flavor='seurat',batch_key='patient');adata.var.to_csv(C/'gene_HVG_audit.tsv',sep='\t');keep=adata.var.highly_variable.values&~adata.var.symbol.isin(heldout).values&~adata.var.symbol.str.startswith(('MT-','RPS','RPL')).values;assert not set(adata.var.loc[keep,'symbol'])&set(heldout)
z=adata[:,keep].copy();print('HVG',z.shape,flush=True);sc.pp.scale(z,max_value=10);sc.tl.pca(z,n_comps=30,svd_solver='arpack',random_state=23);sc.pp.neighbors(z,n_neighbors=30,n_pcs=30,random_state=23);sc.tl.leiden(z,resolution=.5,random_state=23,key_added='cluster');obs['cluster']=z.obs.cluster.values;obs.to_csv(C/'cell_clusters.tsv.gz',sep='\t',index=True);K=obs.cluster.nunique();print('CLUSTERS',K,flush=True)
ct=pd.crosstab(obs.patient,obs.cluster);ct.to_csv(C/'patient_cluster_counts.tsv',sep='\t');summary=[]
for k in ct.columns:
 vv=ct[k];summary.append({'cluster':k,'cells':int(vv.sum()),'patients_ge20':int((vv>=20).sum()),'patients_any':int((vv>0).sum()),'max_patient_fraction':float(vv.max()/vv.sum())})
qc=pd.DataFrame(summary);qc.to_csv(C/'cluster_QC.tsv',sep='\t',index=False)
markers=['C1QA','C1QB','C1QC','CSF1R','LST1','APOE','TREM2','SPP1','FOLR2','LYVE1','MRC1','IL1B','FCN1','S100A8','S100A9','CD1C','FCER1A','CLEC10A','MKI67','TOP2A','ISG15','IFIT3','CD74','HLA-DRA'];wanted=sorted(set(heldout+markers));idx={g:np.flatnonzero(symbols==g) for g in wanted};assert all(len(idx[g])==1 for g in heldout)
rows=[]
# obs position alignment kept; group positions generated explicitly.
obs=obs.reset_index(names='cell_id')
for (pt,k),ix in obs.groupby(['patient','cluster'],observed=True).groups.items():
 ix=np.asarray(list(ix));pb=np.asarray(counts[ix].sum(axis=0)).ravel();lib=pb.sum();n=len(ix)
 for g in wanted:
  if len(idx[g])!=1:continue
  val=float(pb[idx[g][0]]);rows.append({'patient':pt,'cluster':str(k),'sample_type':obs.loc[ix[0],'sample_type'],'n_cells':n,'median_umi':float(np.median(obs.loc[ix,'umi'])),'gene':g,'counts':val,'CPM':val/max(lib,1)*1e6,'log2CPM':float(np.log2(val/max(lib,1)*1e6+1))})
r=pd.DataFrame(rows);r.to_csv(C/'subcluster_pseudobulk_markers.tsv.gz',sep='\t',index=False);r[r.n_cells>=20].groupby(['cluster','gene']).agg(patients=('patient','nunique'),median_log2CPM=('log2CPM','median')).reset_index().to_csv(C/'cluster_marker_patient_medians.tsv',sep='\t',index=False)
prog=r[r.gene.isin(plan['programs']['phagolysosomal'])].groupby(['patient','cluster']).log2CPM.mean().rename('score');tests=[]
for k in sorted(r.cluster.unique()):
 for g in ['CDH23','COPZ2']:
  zz=r[(r.cluster==k)&(r.gene==g)&(r.n_cells>=20)].merge(prog,on=['patient','cluster']);row={'cluster':k,'gene':g,'n':len(zz),'rho':np.nan,'P':np.nan,'partial_rho':np.nan,'P_partial':np.nan}
  if len(zz)>=10 and zz.log2CPM.nunique()>2:
   co=stats.spearmanr(zz.log2CPM,zz.score);row.update(rho=co.statistic,P=co.pvalue);cov=pd.get_dummies(zz.sample_type,drop_first=True,dtype=float).values;A=np.column_stack([np.ones(len(zz)),cov,stats.rankdata(np.log1p(zz.n_cells)),stats.rankdata(np.log1p(zz.median_umi))]);xx=stats.rankdata(zz.log2CPM);yy=stats.rankdata(zz.score);xx-=A@np.linalg.lstsq(A,xx,rcond=None)[0];yy-=A@np.linalg.lstsq(A,yy,rcond=None)[0];df=len(zz)-np.linalg.matrix_rank(A)-1
   if min(np.std(xx),np.std(yy))>1e-9 and df>=3:
    pr=np.corrcoef(xx,yy)[0,1];row.update(partial_rho=pr,P_partial=2*stats.t.sf(abs(pr)*np.sqrt(df/max(1-pr*pr,1e-15)),df))
  tests.append(row)
t=pd.DataFrame(tests);assert len(t)==2*K
for pc,qcname in [('P','FDR_2K'),('P_partial','FDR_partial_2K')]:t[qcname]=stats.false_discovery_control(t[pc].fillna(1).values,method='bh');t.loc[t[pc].isna(),qcname]=np.nan
t.to_csv(C/'within_cluster_associations.tsv',sep='\t',index=False)
# Composition PCs determined without candidate expression.
f=(ct+.5).div((ct+.5).sum(axis=1),axis=0);clr=np.log(f).sub(np.log(f).mean(axis=1),axis=0);pc=PCA(n_components=min(5,K-1),random_state=23).fit(clr);npca=min(len(pc.explained_variance_ratio_),int(np.searchsorted(np.cumsum(pc.explained_variance_ratio_),.8)+1));scores=pd.DataFrame(pc.transform(clr)[:,:npca],index=ct.index,columns=['composition_PC'+str(i+1) for i in range(npca)]);scores.to_csv(C/'composition_PCs.tsv',sep='\t');aa=pd.read_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t');pp=pd.read_csv(O/'external_program_scores.tsv',sep='\t');out=[]
for g in plan['genes']:
 for pname in plan['programs']:
  zz=aa[(aa.celltype=='Macrophage_enriched')&(aa.gene==g)&(aa.n_cells>=20)][['patient','sample_type','n_cells','median_umi','log2CPM']].merge(pp[(pp.celltype=='Macrophage_enriched')&(pp.program==pname)][['patient','score']],on='patient',validate='one_to_one').merge(scores,on='patient');A=np.column_stack([np.ones(len(zz)),pd.get_dummies(zz.sample_type,drop_first=True,dtype=float).values,stats.rankdata(np.log1p(zz.n_cells)),stats.rankdata(np.log1p(zz.median_umi)),zz[scores.columns].values]);xx=stats.rankdata(zz.log2CPM);yy=stats.rankdata(zz.score);xx-=A@np.linalg.lstsq(A,xx,rcond=None)[0];yy-=A@np.linalg.lstsq(A,yy,rcond=None)[0];df=len(zz)-np.linalg.matrix_rank(A)-1;rr=np.corrcoef(xx,yy)[0,1];pv=2*stats.t.sf(abs(rr)*np.sqrt(df/max(1-rr*rr,1e-15)),df);out.append({'gene':g,'program':pname,'n':len(zz),'partial_rho':rr,'P':pv,'df':df,'composition_PCs':npca})
sen=pd.DataFrame(out);sen['FDR_18']=stats.false_discovery_control(sen.P.values,method='bh');sen.to_csv(C/'composition_adjusted_associations.tsv',sep='\t',index=False);(C/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'cells':int(counts.shape[0]),'patients':len(ct),'clusters':int(K),'HVG_after_exclusion':int(keep.sum()),'within_cluster_family':2*K,'composition_PCs':npca,'composition_variance_captured':float(sum(pc.explained_variance_ratio_[:npca])),'within_cluster_significant':int((t.FDR_2K<.05).sum()),'within_cluster_partial_significant':int((t.FDR_partial_2K<.05).sum())},indent=2));print(qc.to_string(index=False));print(t.to_string(index=False));print(sen[sen.program=='phagolysosomal'].to_string(index=False))
