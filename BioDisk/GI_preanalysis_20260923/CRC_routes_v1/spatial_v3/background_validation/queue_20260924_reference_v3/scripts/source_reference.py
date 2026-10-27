from common import *
import os
os.environ['NUMBA_NUM_THREADS']='2';sys.path.insert(0,str(DEEP/'vendor'));import anndata as ad,infercnvpy as cnv
D=Q/'V02_SOURCE';D.mkdir(exist_ok=True);(D/'cells').mkdir(exist_ok=True);(D/'profiles').mkdir(exist_ok=True);samples=json.loads((I/'SAMPLES.json').read_text());co=pd.read_csv(DEEP/'results/cnv/genomic_gene_order.tsv',sep='\t',index_col=0);co.index=co.index.astype(str).str.replace('.','-',regex=False);co=co[~co.index.duplicated(keep=False)];rng=np.random.default_rng(24);meta=[];sources=[];qc=[]
for sample in samples:
 name=sample['sample'];rr=role(sample);X,sy,bc,lib=load_sample(sample);n=len(bc);assert n>0
 panels=C['lineage_markers'];score=[];hits=[]
 for label,gs in panels.items():
  ix=np.flatnonzero(np.isin(sy,gs));z=X[ix].toarray();score.append(np.log1p(z/np.maximum(lib,1)*1e4).mean(axis=0));hits.append((z>0).sum(axis=0))
 score=np.array(score).T;hits=np.array(hits).T;best=score.argmax(axis=1);margin=np.sort(score,axis=1)[:,-1]-np.sort(score,axis=1)[:,-2];labels=np.array(list(panels),dtype=object)[best];valid=(hits[np.arange(n),best]>=2)&(margin>=C['lineage_score_margin'])
 epcore=np.asarray(X[np.isin(sy,['EPCAM','KRT19','KRT20','CDX2'])].sum(axis=0)).ravel()>0;valid&=(labels!='epithelial')|epcore;labels[~valid]='unassigned';labels[labels=='epithelial']='epithelial_unresolved'
 if rr=='liver_normal':labels[labels=='epithelial_unresolved']='biliary_epithelial_like'
 if rr=='colon_normal':labels[labels=='epithelial_unresolved']='normal_colon_epithelial'
 cnvscore=np.full(n,np.nan);pos95=np.zeros(n,bool);pos99=np.zeros(n,bool);refs=[];idx=[];obslabels=[];cal=[]
 epi=np.flatnonzero(labels=='epithelial_unresolved')
 if len(epi)>C['CNV']['epithelial_cap_per_sample']:epi=np.sort(rng.choice(epi,C['CNV']['epithelial_cap_per_sample'],replace=False))
 for label in ['fibroblast','endothelial']:
  v=np.flatnonzero(labels==label)
  if len(v)>=C['CNV']['reference_min_cells']:
   if len(v)>1000:v=rng.choice(v,1000,replace=False)
   v=rng.permutation(v);cut=len(v)//2;idx.extend(v.tolist());obslabels.extend([label+'_train']*cut+['nonreference']*(len(v)-cut));cal.extend([False]*cut+[True]*(len(v)-cut));refs.append(label+'_train')
 before=len(idx);idx.extend(epi.tolist());obslabels.extend(['nonreference']*len(epi));cal.extend([False]*len(epi));idx=np.array(idx,dtype=int);cnv_status='insufficient_reference_or_epithelial';cut95=np.nan;cut99=np.nan
 if len(epi)>=20 and refs:
  unique=~pd.Index(sy).duplicated(keep=False);vi=np.flatnonzero(unique&np.isin(sy,co.index));mean=np.asarray(X[vi][:,idx].mean(axis=1)).ravel();vi=vi[mean>=.1]
  if len(vi)>=C['CNV']['min_genes']:
   var=co.reindex(sy[vi]);obs=pd.DataFrame({'ref':obslabels,'cal':cal,'broad':labels[idx]},index=[name+'~'+b for b in bc[idx]]);Y=X[vi][:,idx].T.tocsr().multiply((1e4/lib[idx])[:,None]).tocsr();Y.data=np.log1p(Y.data);aa=ad.AnnData(Y,obs=obs,var=var);cnv.tl.infercnv(aa,reference_key='ref',reference_cat=refs,window_size=100,step=10,dynamic_threshold=1.5,chunksize=500,n_jobs=1);z=np.asarray(abs(aa.obsm['X_cnv']).mean(axis=1)).ravel();cnvscore[idx]=z
   cut95=max(np.quantile(z[(obs.broad==ct).to_numpy()&obs.cal.to_numpy()],.95) for ct in ['fibroblast','endothelial'] if ((obs.broad==ct)&obs.cal).any());cut99=max(np.quantile(z[(obs.broad==ct).to_numpy()&obs.cal.to_numpy()],.99) for ct in ['fibroblast','endothelial'] if ((obs.broad==ct)&obs.cal).any());pos95[epi]=cnvscore[epi]>cut95;pos99[epi]=cnvscore[epi]>cut99;labels[pos95]='epithelial_CNV_supported';cnv_status='RNA_CNV_calibrated'
  else:cnv_status='insufficient_CNV_genes'
 cell=pd.DataFrame({'barcode':bc,'label':labels,'lineage_margin':margin,'RNA_CNV_score':cnvscore,'CNV95':pos95,'CNV99':pos99,'UMI':lib});cell.to_csv(D/'cells'/(name+'.tsv.gz'),sep='\t',index=False)
 names=[];profs=[];measure=~pd.Index(sy).duplicated(keep=False)
 for label in sorted(set(labels)-{'unassigned'}):
  k=np.flatnonzero(labels==label)
  for gene in C['targets']:
   g=np.flatnonzero(sy==gene)
   if len(g)==1:sources.append({'sample':name,'dataset':sample['dataset'],'donor_label':sample['patient_label_provisional'],'role':rr,'label':label,'cells':len(k),'gene':gene,'mean_CP10K':float(np.asarray(X[g][:,k].multiply(1e4/lib[k]).mean(axis=1)).ravel()[0]),'detected_pct':float(np.mean(X[g][:,k].toarray().ravel()>0)*100)})
  if len(k)<20:continue
  profile=np.asarray(X[:,k].multiply(1e6/lib[k]).mean(axis=1)).ravel();names.append(label);profs.append(profile[measure]);meta.append({'sample':name,'dataset':sample['dataset'],'donor_label':sample['patient_label_provisional'],'role':rr,'label':label,'cells':len(k),'file':name+'.npz'})
 if profs:np.savez_compressed(D/'profiles'/(name+'.npz'),genes=sy[measure].astype(str),labels=np.array(names,dtype=str),profiles=np.stack(profs,axis=1))
 qc.append({'sample':name,'role':rr,'QC_cells':n,'confident_lineage_cells':int((labels!='unassigned').sum()),'CNV_status':cnv_status,'CNV95':int(pos95.sum()),'CNV99':int(pos99.sum()),'cutoff95':cut95,'cutoff99':cut99});pd.DataFrame(qc).to_csv(D/'sample_QC.tsv',sep='\t',index=False);print(name,cnv_status,int(pos95.sum()),flush=True)
pd.DataFrame(meta).to_csv(D/'profile_metadata.tsv',sep='\t',index=False);pd.DataFrame(sources).to_csv(D/'source_expression.tsv',sep='\t',index=False);m=pd.DataFrame(meta);keep=m[(m.cells>=50)&m.role.str.startswith('liver')];required=C['reference_gate']['required'];gate={ct:int(keep[keep.label==ct].donor_label.nunique()) for ct in required};allowed=all(v>=3 for v in gate.values());dump(D/'GATE.json',{'allowed':allowed,'donor_labels_per_required_type':gate,'reason':'minimum3 donor labels with50 cells per profile' if allowed else 'reference coverage insufficient','scope':'Exploratory reference; donor independence/clinical site unresolved; RNA-CNV support is not DNA confirmation'})
if allowed:
 allgenes=None;cache={}
 for name in keep['sample'].unique():
  z=np.load(D/'profiles'/(name+'.npz'));cache[name]=z;allgenes=set(z['genes']) if allgenes is None else allgenes&set(z['genes'])
 genes=np.array(sorted(allgenes));types=required+['epithelial_unresolved','biliary_epithelial_like'];profiles=[];used=[]
 for ct in types:
  rows=keep[keep.label==ct];donors=[]
  for donor,v in rows.groupby('donor_label'):
   vals=[]
   for _,r in v.iterrows():
    z=cache[r['sample']];gi=pd.Index(z['genes']).get_indexer(genes);ci=np.flatnonzero(z['labels']==ct)[0];vals.append(z['profiles'][gi,ci])
   donors.append(np.mean(vals,axis=0))
  if len(donors)>=3:used.append(ct);profiles.append(np.mean(donors,axis=0))
 np.savez_compressed(D/'reference.npz',genes=genes.astype(str),types=np.array(used,dtype=str),CPM=np.stack(profiles,axis=1))
dump(D/'RESULT.json',{'status':'computed_pending_review','samples':len(qc),'mixture_gate_allowed':allowed,'malignancy_scope':'RNA-inferred CNV only','clinical_inference':'not_ready'})
