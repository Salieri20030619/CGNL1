from pathlib import Path
import os,json,gzip,datetime,hashlib,sys,gc,importlib.metadata
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse,io
from sklearn.metrics import adjusted_rand_score,balanced_accuracy_score,confusion_matrix
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];B=R.parent;C=json.loads((R/'SPECIFICATION.json').read_text());EX=set(C['excluded_target_panel']);EX10=EX|set(C['extra_exclusion_sensitivity']);SEED=20261002
sc.settings.n_jobs=4
plt.rcParams.update({'svg.fonttype':'none','font.size':9})
def lognorm(X):
 X=X.tocsr().astype(np.float32);lib=np.asarray(X.sum(axis=1)).ravel();assert (lib>0).all();X=sparse.diags(1e4/lib)@X;X=X.tocsr();np.log1p(X.data,out=X.data);return X

def class_profiles(a,exclude_donor=None):
 rows=[];keys=[]
 m=a.obs
 for (donor,label),g in m.groupby(['donor','label'],observed=True):
  if exclude_donor is not None and donor==exclude_donor:continue
  ii=m.index.get_indexer(g.index);rows.append(np.asarray(a.X[ii].mean(axis=0)).ravel());keys.append((donor,label))
 assert len(rows)>0
 rows=np.stack(rows);key=pd.DataFrame(keys,columns=['donor','label']);labels=sorted(key.label.unique());P=np.stack([rows[key.label.eq(l)].mean(axis=0) for l in labels]);return labels,P,key

def choose_features(a,exclude_donor=None,allowed=None):
 labels,P,key=class_profiles(a,exclude_donor);genes=a.var_names.to_numpy(dtype=str);valid=np.array([g not in EX10 and not g.startswith(('MT-','RPL','RPS')) for g in genes])
 if allowed is not None:valid &= np.isin(genes,list(allowed))
 selected=set()
 for j,label in enumerate(labels):
  diff=P[j]-(P.sum(axis=0)-P[j])/(len(labels)-1);cand=np.flatnonzero(valid&(P[j]>.05));selected.update(cand[np.argsort(diff[cand])[-150:]].tolist())
 ii=np.array(sorted(selected));assert len(ii)>=300
 return labels,P[:,ii],genes[ii],ii

def correlations(Q,P):
 Q=np.asarray(Q,dtype=np.float64);P=np.asarray(P,dtype=np.float64);Q-=Q.mean(axis=1,keepdims=True);P=P-P.mean(axis=1,keepdims=True)
 return Q@P.T/np.maximum(np.linalg.norm(Q,axis=1)[:,None]*np.linalg.norm(P,axis=1)[None,:],1e-12)

def predict_matrix(X,centroids,labels):
 results=[];best=[];margins=[]
 for start in range(0,X.shape[0],2048):
  block=X[start:start+2048].toarray() if sparse.issparse(X) else X[start:start+2048];s=correlations(block,centroids);arg=s.argmax(axis=1);ss=np.sort(s,axis=1);results.extend(np.array(labels)[arg]);best.extend(ss[:,-1]);margins.extend(ss[:,-1]-ss[:,-2])
 return np.array(results),np.array(best),np.array(margins)

def prepare_refs():
 refs={}
 # Disease-matched author reference; expression is adjusted rather than raw UMI.
 f=B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad';a=ad.read_h5ad(f);m=a.obs.copy();mapping={'Epithelial':'epithelial','Myeloid':'myeloid','T-NK-ILC':'T_NK','B':'B','Plasma':'plasma','Endothelial':'endothelial','Hepatocyte':'hepatocyte','Stromal':'fibroblast'};m['label']=m.Cell_type.astype(str).map(mapping);m.loc[m.Cell_subtype.astype(str).eq('Cholangiocyte'),'label']='cholangiocyte';m.loc[m.Cell_subtype.astype(str).isin(['Pericyte','SMC']),'label']='mural';m['donor']=m.Patient.astype(str);assert m.label.notna().all()
 keep=~a.var_names.isin(EX);X=lognorm(a.X[:,keep]);refs['multiome']=ad.AnnData(X,obs=m[['label','donor']],var=pd.DataFrame(index=a.var_names[keep]));del a,X
 # Independent normal-liver reference, known log2 normalized author expression.
 p=B/'spatial_v3/background_validation/reference_supplement_20260924';m=pd.read_csv(p/'GSE115469_CellClusterType.txt.gz',sep='\t').set_index('CellName');chunks=[];genes=[];cols=None
 for t in pd.read_csv(p/'GSE115469_Data.csv.gz',index_col=0,chunksize=256):
  if cols is None:cols=t.columns;m=m.reindex(cols);assert m.CellType.notna().all()
  t=t.loc[~t.index.isin(EX)];z=t.to_numpy(np.float32);assert np.isfinite(z).all() and z.min()>=0 and z.max()<40;chunks.append(sparse.csr_matrix(np.expm1(z*np.log(2))));genes.extend(t.index.astype(str))
 X=sparse.vstack(chunks,format='csr').T.tocsr();del chunks
 def broad(s):
  if 'Hepatocyte' in s:return 'hepatocyte'
  if 'LSEC' in s or 'endothelial' in s:return 'endothelial'
  if 'Macrophage' in s:return 'myeloid'
  if 'T_Cells' in s or 'NK-like' in s:return 'T_NK'
  if 'Plasma' in s:return 'plasma'
  if 'B_Cells' in s:return 'B'
  if 'Cholangio' in s:return 'cholangiocyte'
  if 'Stellate' in s:return 'fibroblast'
  if 'Erythroid' in s:return 'erythroid'
  raise ValueError(s)
 m['label']=m.CellType.map(broad);m['donor']=m.Sample.astype(str);refs['normal_liver']=ad.AnnData(lognorm(X),obs=m[['label','donor']],var=pd.DataFrame(index=genes));del X
 for name,a in refs.items():
  assert a.var_names.is_unique and a.obs_names.is_unique and not set(a.var_names)&EX
  a.write_h5ad(R/'references'/(name+'_target_excluded.h5ad'),compression='gzip');validation=[];preds=[]
  for donor in sorted(a.obs.donor.unique()):
   labels,P,gs,ii=choose_features(a,donor);mask=a.obs.donor.eq(donor).to_numpy();pred,score,margin=predict_matrix(a.X[mask][:,ii],P,labels);truth=a.obs.loc[mask,'label'].to_numpy();validation.append(dict(reference=name,donor=donor,n=int(mask.sum()),accuracy=float(np.mean(pred==truth)),balanced_accuracy=float(balanced_accuracy_score(truth,pred)),endo_cells=int((truth=='endothelial').sum()),endo_recall=float(np.mean(pred[truth=='endothelial']=='endothelial')) if (truth=='endothelial').any() else np.nan))
   preds.extend([dict(donor=donor,truth=t,predicted=p) for t,p in zip(truth,pred)])
  pd.DataFrame(validation).to_csv(R/'references'/(name+'_LODO.tsv'),sep='\t',index=False);pd.crosstab(pd.DataFrame(preds).truth,pd.DataFrame(preds).predicted).to_csv(R/'references'/(name+'_LODO_confusion.tsv'),sep='\t');print('REFERENCE_READY',name,a.shape,flush=True)
 return refs

def load_cohort(ds):
 I=B/'spatial_v3/background_validation/run_20260924T052010Z';jobs=[]
 if ds=='GSE285991':
  for s in json.loads((I/'SAMPLES.json').read_text()):
   if s['dataset']==ds:jobs.append((s['sample'],s['patient_label_provisional'],[I/x for x in s['files']]))
 else:
  for s in json.loads((B/'cohort_verification_20260924/GSE245552_LM_SAMPLES.json').read_text()):jobs.append((s['sample'],ds+':'+s['donor'],[B/'cohort_verification_20260924/GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]))
 objects=[];qc=[]
 for sample,donor,paths in jobs:
  fs={k:next(p for p in paths if k in p.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(fs['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(fs['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
  with gzip.open(fs['matrix'],'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
  X.eliminate_zeros();assert X.shape==(len(sy),len(bc));lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;good=(lib>=200)&(ng>=100)&(mt<=35)
  keep=~np.isin(sy,list(EX));X=X[keep][:,good];sy=sy[keep];unique,inverse=np.unique(sy,return_inverse=True);collapse=sparse.csr_matrix((np.ones(len(sy)),(inverse,np.arange(len(sy)))),shape=(len(unique),len(sy)));X=(collapse@X).T.tocsr().astype(np.float32);meta=pd.DataFrame({'sample':sample,'donor':donor,'barcode':bc[good],'UMI_QC':lib[good],'genes_QC':ng[good],'mito_QC':mt[good]},index=[sample+':'+s for s in bc[good]]);objects.append(ad.AnnData(X,obs=meta,var=pd.DataFrame(index=unique)));qc.append(dict(sample=sample,donor=donor,input_cells=len(bc),QC_cells=int(good.sum())));print('LOADED',ds,sample,int(good.sum()),flush=True)
 a=ad.concat(objects,join='outer',merge='same');del objects;a.X=a.X.tocsr();assert not set(a.var_names)&EX;a.X=lognorm(a.X);return a,pd.DataFrame(qc)

def cluster(a,extra=False,umap=False):
 keep=~a.var_names.isin(EX10 if extra else EX);q=a[:,keep].copy();sc.pp.highly_variable_genes(q,n_top_genes=3000,flavor='seurat',batch_key='sample');q=q[:,q.var.highly_variable].copy();assert not set(q.var_names)&(EX10 if extra else EX);sc.pp.scale(q,max_value=10);sc.tl.pca(q,n_comps=30,svd_solver='arpack',random_state=SEED);sc.pp.neighbors(q,n_neighbors=20,n_pcs=30,metric='cosine',random_state=SEED)
 for res,key in [(.5,'cluster05'),(1.0,'cluster10')]:sc.tl.leiden(q,resolution=res,key_added=key,random_state=SEED,n_iterations=2)
 if umap:sc.tl.umap(q,min_dist=.3,maxiter=200,random_state=SEED)
 return q

def broadlabel(x):return 'immune' if x in ['myeloid','T_NK','B','plasma','erythroid'] else ('stroma' if x in ['fibroblast','mural'] else x)
def categorical_plot(ax,xy,values,title,order=None):
 values=np.asarray(values).astype(str);names=order or sorted(set(values));colors=plt.cm.tab20(np.linspace(0,1,max(len(names),2)))
 for i,name in enumerate(names):
  k=values==name;ax.scatter(xy[k,0],xy[k,1],s=1,c=[colors[i]],label=name,rasterized=True,alpha=.65)
 ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.legend(markerscale=5,fontsize=6,loc='best',ncol=2)

def process(ds,refs):
 D=R/'results'/ds;D.mkdir();a,qc=load_cohort(ds);qc.to_csv(D/'sample_QC.tsv',sep='\t',index=False);print('CLUSTER_START',ds,a.shape,flush=True);q=cluster(a,umap=True);pd.Series(q.var_names,name='gene').to_csv(D/'primary_HVG.tsv',sep='\t',index=False);print('PRIMARY_CLUSTER_DONE',ds,flush=True)
 alt=cluster(a,extra=True);pd.Series(alt.var_names,name='gene').to_csv(D/'extra_exclusion_HVG.tsv',sep='\t',index=False);out=a.obs.copy();out['cluster05']=q.obs.cluster05.to_numpy();out['cluster10']=q.obs.cluster10.to_numpy();out['extra_cluster05']=alt.obs.cluster05.to_numpy();out['extra_cluster10']=alt.obs.cluster10.to_numpy();ari=adjusted_rand_score(out.cluster05,out.extra_cluster05);del alt;gc.collect()
 for name,ref in refs.items():
  labels,P,gs,ii=choose_features(ref,allowed=set(a.var_names));assert not set(gs)&EX10;pd.DataFrame(P,index=labels,columns=gs).to_csv(D/(name+'_centroids.tsv'),sep='\t');pd.Series(gs,name='gene').to_csv(D/(name+'_features.tsv'),sep='\t',index=False)
  ids=a.var_names.get_indexer(gs);pred,score,margin=predict_matrix(a.X[:,ids],P,labels);out[name+'_prediction']=pred;out[name+'_score']=score;out[name+'_margin']=margin
  cm=[]
  for cluster_id,g in out.groupby('cluster05',observed=True):
   ids2=out.index.get_indexer(g.index);mean=np.asarray(a.X[ids2][:,ids].mean(axis=0)).ravel();ss=correlations(mean[None,:],P)[0];cm.append(dict(cluster=cluster_id,prediction=labels[int(ss.argmax())],score=float(ss.max()),margin=float(np.sort(ss)[-1]-np.sort(ss)[-2])))
  pd.DataFrame(cm).to_csv(D/(name+'_cluster_reference.tsv'),sep='\t',index=False)
 out['multiome_broad']=out.multiome_prediction.map(broadlabel);out['normal_broad']=out.normal_liver_prediction.map(broadlabel);out['reference_agreement']=np.where(out.multiome_broad==out.normal_broad,'both_'+out.multiome_broad,'discordant')
 # Only now join the old identity flags; no target values enter this table.
 out['old_group']='other_QC_cells';cache=json.loads((B/'threshold_ladder_20261002_v1/CACHE_INDEX.json').read_text())
 for rec in cache:
  if rec['dataset']!=ds:continue
  with np.load(rec['path']) as z:
   ids=pd.Index([rec['sample']+':'+b for b in z['barcode']]);present=ids.isin(out.index);qcm=(z['qlib']>=500)&(z['qng']>=200)&(z['qmt']<=25);group=np.where(qcm,np.where(z['hit']==1,'one_marker_removed','two_plus_retained'),'broad_only_QC');out.loc[ids[present],'old_group']=group[present]
 out['UMAP1']=q.obsm['X_umap'][:,0];out['UMAP2']=q.obsm['X_umap'][:,1];out.to_csv(D/'cell_identity_evidence.tsv.gz',sep='\t',index=True,index_label='cell_id')
 # Cluster marker differences are descriptive, not tests using cells as replicates.
 markers=[];features=a.var_names.to_numpy(dtype=str);globalmean=np.asarray(a.X.mean(axis=0)).ravel();globalhits=a.X.getnnz(axis=0)
 for cl,g in out.groupby('cluster05',observed=True):
  k=out.index.get_indexer(g.index);Z=a.X[k];mean=np.asarray(Z.mean(axis=0)).ravel();pct=Z.getnnz(axis=0)/len(k);other=(globalmean*a.n_obs-mean*len(k))/max(a.n_obs-len(k),1);delta=mean-other;valid=(pct>=.1)&~np.char.startswith(features,'MT-')&~np.char.startswith(features,'RPL')&~np.char.startswith(features,'RPS');ix=np.flatnonzero(valid);ix=ix[np.argsort(delta[ix])[-30:][::-1]]
  for j in ix:markers.append(dict(cluster=cl,gene=features[j],mean_log=mean[j],mean_difference=delta[j],fraction_detected=pct[j],cells=len(k)))
 pd.DataFrame(markers).to_csv(D/'cluster_markers_descriptive.tsv',sep='\t',index=False)
 pd.crosstab(out.old_group,out.reference_agreement).to_csv(D/'old_group_reference_agreement.tsv',sep='\t');pd.crosstab(out.old_group,out.cluster05).to_csv(D/'old_group_cluster_counts.tsv',sep='\t');pd.crosstab(out.cluster05,out.donor).to_csv(D/'cluster_donor_counts.tsv',sep='\t');pd.crosstab(out.cluster05,out.multiome_prediction,normalize='index').to_csv(D/'cluster_reference_fractions.tsv',sep='\t')
 per=[]
 for (donor,group),g in out[out.old_group!='other_QC_cells'].groupby(['donor','old_group'],observed=True):
  per.append(dict(donor=donor,group=group,cells=len(g),both_immune=float(g.reference_agreement.eq('both_immune').mean()),both_endo=float(g.reference_agreement.eq('both_endothelial').mean()),discordant=float(g.reference_agreement.eq('discordant').mean()),multiome_endo=float(g.multiome_prediction.eq('endothelial').mean()),normal_endo=float(g.normal_liver_prediction.eq('endothelial').mean())))
 pd.DataFrame(per).to_csv(D/'patient_old_group_reference_fractions.tsv',sep='\t',index=False)
 fig,ax=plt.subplots(2,2,figsize=(14,11));xy=q.obsm['X_umap'];categorical_plot(ax[0,0],xy,out.cluster05,'Target-excluded clustering');categorical_plot(ax[0,1],xy,out.multiome_prediction,'CRCLM author reference');categorical_plot(ax[1,0],xy,out.normal_liver_prediction,'Normal-liver author reference');categorical_plot(ax[1,1],xy,out.old_group,'Previous candidate groups (overlay only)');fig.suptitle(ds+' | identity review without CGNL1 / junction panel',fontsize=14);fig.tight_layout(rect=[0,0,1,.97]);fig.savefig(D/'identity_umap.png',dpi=160);plt.close(fig)
 embedding=ad.AnnData(obs=out.copy());embedding.obsm['X_umap']=q.obsm['X_umap'];embedding.obsm['X_pca']=q.obsm['X_pca'];embedding.write_h5ad(D/'identity_embeddings.h5ad',compression='gzip')
 summary={'exit_code':0,'dataset':ds,'QC_cells':a.n_obs,'patients':out.donor.nunique(),'clusters05':out.cluster05.nunique(),'clusters10':out.cluster10.nunique(),'extra_marker_exclusion_ARI05':ari,'target_and_panel_features_used':0,'reference_original_identity_features_used':0,'association_tests_run':False,'old_group_reference_counts':pd.crosstab(out.old_group,out.reference_agreement).to_dict(),'finished':datetime.datetime.now().astimezone().isoformat()};(D/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print('COHORT_COMPLETE',ds,json.dumps(summary),flush=True);del a,q,out;gc.collect()
refs=prepare_refs()
for ds in ['GSE245552','GSE285991']:process(ds,refs)
(R/'ENVIRONMENT.json').write_text(json.dumps({k:importlib.metadata.version(k) for k in ['scanpy','anndata','numpy','scipy','pandas','leidenalg','igraph']},indent=2));(R/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'cohorts':2,'references':2,'scope':'target-excluded clustering and reference identity trial only','association_tests_run':False,'finished':datetime.datetime.now().astimezone().isoformat()},indent=2));print('ALL_COMPLETE',flush=True)
