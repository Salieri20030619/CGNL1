
from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
from sklearn.metrics import balanced_accuracy_score
from threadpoolctl import threadpool_limits
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'A04_subtype_check_20261002_v1';D=R/'identity';D.mkdir(parents=True,exist_ok=True)
EX={'CGNL1','CLDN5','ESAM','TJP1','OCLN','PECAM1','VWF','KDR','EMCN','CDH5'}
panels={'sinusoidal':['CLEC4G','CLEC4M','STAB2','FCGR2B','CLEC1B','FCN2','FCN3','OIT3'],'arterial':['GJA5','EFNB2','SOX17','SEMA3G','DLL4'],'venous':['ACKR1','SELP','NR2F2'],'capillary':['CA4','RGCC','RBP7','PLVAP'],'lymphatic':['PROX1','PDPN','CCL21','FLT4'],'tip_state':['ESM1','ANGPT2','APLN','CXCR4'],'cycling':['MKI67','TOP2A','UBE2C'],'immune':['PTPRC','LST1','TYROBP','CD3D']}
spec={'created':datetime.datetime.now().astimezone().isoformat(),'scope':'A04 diagnosis in two reviewed cohorts; does not claim whole-plan completion or independent validation','excluded_genes':sorted(EX),'marker_panels':panels,'references':['GSE115469 author 3 endothelial labels','CRCLM multiome author 4 endothelial subtypes'],'annotation':'patient-balanced nearest centroid Pearson; top100 reference-only contrasts per class; no query target outcome used; no confidence-based exclusion in broad primary','coarse_map':{'SEC':'sinusoidal_like','Lymphatic endo':'lymphatic_like','Tip-like endo':'other_blood_endo_like','Stalk-like endo':'other_blood_endo_like'},'normal_reference_limit':'no lymphatic class; cannot adjudicate absence of lymphatic cells','clustering':'separate cohort; HVG1500 sample batch, PCA20, k15 cosine, Leiden0.3/0.6 seed20261002; excluded target/panel/original identity genes; no integration for outcome inference','planned_tests':'fixed n10 lesion primary18: study, study+sinusoidal fraction, study+lymphatic fraction, study+two fractions, study+two fractions+depth; broad n5 per subtype lesion, sensitivities n3/n10; matched caliper.25 total3 and10 pairs with same coarse subtype restriction, compare depth-only on common eligible patients; primary leave-one-out and raw-detection summaries','FDR':'composition models one family; subtype tests per cell threshold family; paired matching family; influence diagnostic family; all target P/FDR retained','source_urls':['https://www.nature.com/articles/s41467-018-06318-7','https://www.nature.com/articles/s41586-019-1373-2','https://www.nature.com/articles/s41591-024-03376-x'],'status':'exploratory reference-supported classification, not independent blind subtype adjudication'}
(R/'SPECIFICATION.json').write_text(json.dumps(spec,indent=2))
refs={}
a=ad.read_h5ad(B/'identity_review_target_blind_20261002_v1/references/multiome_target_excluded.h5ad');src=ad.read_h5ad(B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad',backed='r');a=a[a.obs.label.eq('endothelial')].copy();a.obs['subtype']=src.obs.loc[a.obs_names,'Cell_subtype'].astype(str);src.file.close();refs['disease']=a
a=ad.read_h5ad(B/'identity_review_target_blind_20261002_v1/references/normal_liver_target_excluded.h5ad');src=pd.read_csv(B/'spatial_v3/background_validation/reference_supplement_20260924/GSE115469_CellClusterType.txt.gz',sep='\t').set_index('CellName');a=a[a.obs.label.eq('endothelial')].copy();a.obs['subtype']=src.loc[a.obs_names,'CellType'].astype(str);refs['normal']=a
def model(a,omit=None,allowed=None):
 mask=np.ones(a.n_obs,bool) if omit is None else (a.obs.donor!=omit).to_numpy();t=a.obs[mask];xx=a.X[mask];rows=[];keys=[]
 for (donor,label),g in t.groupby(['donor','subtype'],observed=True):
  ix=t.index.get_indexer(g.index);rows.append(np.asarray(xx[ix].mean(axis=0)).ravel());keys.append(label)
 labs=sorted(set(keys));P=np.stack([np.mean([x for x,k in zip(rows,keys) if k==l],axis=0) for l in labs]);genes=a.var_names.to_numpy();ok=np.array([g not in EX and not g.startswith(('MT-','RPL','RPS')) and (allowed is None or g in allowed) for g in genes]);sel=set()
 for j in range(len(labs)):
  diff=P[j]-(P.sum(axis=0)-P[j])/(len(labs)-1);ids=np.flatnonzero(ok&(P[j]>.05));sel.update(ids[np.argsort(diff[ids])[-100:]].tolist())
 ids=np.array(sorted(sel));return labs,P[:,ids],genes[ids],ids
def predict(X,P,labs):
 X=X.toarray() if sparse.issparse(X) else np.asarray(X);x=X-X.mean(axis=1,keepdims=True);p=P-P.mean(axis=1,keepdims=True);z=x@p.T/np.maximum(np.linalg.norm(x,axis=1)[:,None]*np.linalg.norm(p,axis=1),1e-12);args=z.argmax(axis=1);ss=np.sort(z,axis=1);return np.array(labs)[args],z,ss[:,-1]-ss[:,-2]
validation=[];conf=[]
for name,a in refs.items():
 a.obs[['donor','subtype']].to_csv(D/(name+'_reference_labels.tsv'),sep='\t')
 for donor in a.obs.donor.unique():
  labs,P,genes,ids=model(a,donor);k=(a.obs.donor==donor).to_numpy();pr,z,ma=predict(a.X[k][:,ids],P,labs);truth=a.obs.loc[k,'subtype'].to_numpy();validation.append(dict(reference=name,donor=donor,n=int(k.sum()),accuracy=float(np.mean(pr==truth)),balanced_accuracy=float(balanced_accuracy_score(truth,pr)),training_labels='|'.join(labs)))
  conf.extend(dict(reference=name,donor=donor,truth=t,predicted=p) for t,p in zip(truth,pr))
pd.DataFrame(validation).to_csv(D/'reference_LODO.tsv',sep='\t',index=False);cf=pd.DataFrame(conf);cf.groupby(['reference','truth','predicted']).size().rename('cells').reset_index().to_csv(D/'reference_LODO_confusion.tsv',sep='\t',index=False)
cache=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text());results=[];markers=[];scores=[];mix=[]
for ds in ['GSE245552','GSE285991']:
 arrays=[]
 for rec in cache:
  if rec['dataset']!=ds or rec['selected_cells']==0:continue
  genes=pd.read_csv(rec['genes'],sep='\t').gene.to_numpy(dtype=str);X=sparse.load_npz(rec['path']).T.tocsr();g=pd.read_csv(rec['cells'],sep='\t',dtype={'cluster05':str}).set_index('cell_id');g['donor']=rec['donor'];assert X.shape[0]==len(g);keep=(~pd.Index(genes).duplicated(keep=False))&(~np.isin(genes,list(EX)));arrays.append(ad.AnnData(X[:,keep],obs=g[['sample','donor','barcode']],var=pd.DataFrame(index=genes[keep])))
 a=ad.concat(arrays,join='inner');assert a.obs_names.is_unique and not set(a.var_names)&EX;sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);out=a.obs.copy()
 for name,ref in refs.items():
  labs,P,genes,ids=model(ref,allowed=set(a.var_names));assert not set(genes)&EX;pd.DataFrame(P,index=labs,columns=genes).to_csv(D/(ds+'_'+name+'_centroids.tsv'),sep='\t');ii=a.var_names.get_indexer(genes);pr,z,margin=predict(a.X[:,ii],P,labs);out[name+'_subtype']=pr;out[name+'_score']=z.max(axis=1);out[name+'_margin']=margin
 out['coarse_subtype']=out.disease_subtype.map(spec['coarse_map']);assert out.coarse_subtype.notna().all()
 for key,gs in panels.items():
  present=[g for g in gs if g in a.var_names];out[key+'_score']=np.asarray(a[:,present].X.mean(axis=1)).ravel();out[key+'_hits']=np.asarray((a[:,present].X>0).sum(axis=1)).ravel()
 # Label certainty is descriptive only, never an exclusion rule in broad branch.
 out['sinusoid_reference_agreement']=(out.coarse_subtype=='sinusoidal_like')&out.normal_subtype.isin(['Periportal_LSECs','Central_venous_LSECs'])
 q=a.copy();sc.pp.highly_variable_genes(q,n_top_genes=min(1500,q.n_vars),batch_key='sample',flavor='seurat');q=q[:,q.var.highly_variable].copy();sc.pp.scale(q,max_value=10);sc.tl.pca(q,n_comps=20,svd_solver='arpack',random_state=20261002);sc.pp.neighbors(q,n_neighbors=15,n_pcs=20,metric='cosine',random_state=20261002);sc.tl.leiden(q,resolution=.3,key_added='cluster03',random_state=20261002);sc.tl.leiden(q,resolution=.6,key_added='cluster06',random_state=20261002);out['cluster03']=q.obs.cluster03;out['cluster06']=q.obs.cluster06
 for cluster,g in out.groupby('cluster03',observed=True):
  ix=out.index.get_indexer(g.index);rest=np.flatnonzero(~out.index.isin(g.index));mean=np.asarray(a.X[ix].mean(axis=0)).ravel();other=np.asarray(a.X[rest].mean(axis=0)).ravel();det=np.asarray((a.X[ix]>0).mean(axis=0)).ravel();delta=mean-other;valid=np.array([not s.startswith(('MT-','RPL','RPS')) for s in a.var_names]);top=np.flatnonzero(valid)[np.argsort(delta[valid])[-30:][::-1]]
  for j in top:markers.append(dict(dataset=ds,cluster=cluster,gene=a.var_names[j],mean=mean[j],detection_fraction=det[j],mean_difference=delta[j]))
  mix.append(dict(dataset=ds,cluster=cluster,cells=len(g),donors=g.donor.nunique(),largest_donor_fraction=g.donor.value_counts(normalize=True).max(),**{l:float(g.coarse_subtype.eq(l).mean()) for l in spec['coarse_map'].values()}))
  for key in panels:scores.append(dict(dataset=ds,cluster=cluster,panel=key,score=g[key+'_score'].mean(),mean_detected_markers=g[key+'_hits'].mean()))
 out['dataset']=ds;out.to_csv(D/(ds+'_cell_subtype.tsv.gz'),sep='\t',index_label='cell_id');results.append(out);print('ANNOTATED',ds,len(out),out.coarse_subtype.value_counts().to_dict(),flush=True)
pd.DataFrame(markers).to_csv(D/'cluster_markers_descriptive.tsv',sep='\t',index=False);pd.DataFrame(scores).to_csv(D/'cluster_panel_scores.tsv',sep='\t',index=False);pd.DataFrame(mix).to_csv(D/'cluster_composition.tsv',sep='\t',index=False)
combined=pd.concat(results);combined.groupby(['dataset','donor','coarse_subtype']).size().rename('cells').reset_index().to_csv(D/'donor_subtype_counts.tsv',sep='\t',index=False)
freeze={'time':datetime.datetime.now().astimezone().isoformat(),'target_outcomes_tested':False,'excluded':sorted(EX),'annotation_cells':len(combined),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.glob('*') if p.is_file()}}
(R/'IDENTITY_FREEZE.json').write_text(json.dumps(freeze,indent=2));print('IDENTITY_COMPLETE',flush=True)
