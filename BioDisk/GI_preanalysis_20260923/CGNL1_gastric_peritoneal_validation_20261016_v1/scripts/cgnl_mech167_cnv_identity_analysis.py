from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
import json,hashlib,datetime,numpy as np,pandas as pd,anndata as ad,infercnvpy as cnv
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path.cwd()/'mechanism_investigation_v2';spec=json.loads((O/'CNV_IDENTITY_SPECIFICATION_165.json').read_text(encoding='utf-8'));inputs=json.loads((O/'metadata/CNV_INPUTS_READY_166.json').read_text(encoding='utf-8'));ann=json.loads((O/'metadata/CNV_ANNOTATION_READY_164.json').read_text(encoding='utf-8'));assert hashlib.sha256(Path(ann['gene_order']).read_bytes()).hexdigest()==ann['gene_order_sha256'];genes=pd.read_csv(ann['gene_order'],sep='\t').set_index('gene_id');assert genes.index.is_unique;root=O/'results/CNV_identity';root.mkdir(exist_ok=True);records=[];summaries=[];similarities=[]
group_order=['immune_like','endothelial_like','fibroblast_like','mural_like','overlapping_stromal','unclassified','strict_mesothelial','conservative_epithelial']
def serial(x):
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,np.ndarray):return x.tolist()
 raise TypeError(type(x).__name__)
for item in inputs['cases']:
 p=Path(item['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'];a=ad.read_h5ad(p);mapped=genes.reindex(a.var_names);det=np.asarray((a.layers['counts']>0).sum(0)).ravel();keep=mapped.chromosome.notna().to_numpy()&(det>=20);initial=a.n_vars;a=a[:,keep].copy();mapped=genes.loc[a.var_names];a.var['chromosome']=mapped.chromosome.astype(str);a.var['start']=mapped.start.astype(np.int64);a.var['end']=mapped.end.astype(np.int64);assert not set(a.var.gene_symbol)&set(spec['gene_exclusions']);assert a.n_vars>5000;sample=item['dataset']+'_'+item['sample'];D=root/sample;D.mkdir(exist_ok=True);a.var.to_csv(D/'mapped_genes.tsv',sep='\t',encoding='utf-8');case={'dataset':item['dataset'],'sample':item['sample'],'input_sha256':item['sha256'],'cells':a.n_obs,'input_genes':initial,'mapped_expressed_autosomal_genes':a.n_vars,'models':[]};profiles={}
 for model,refs in [('primary',spec['primary_reference']),('immune_reference',spec['sensitivity_reference'])]:
  b=a.copy();key='cnv_'+model;print('CNV_RUN',sample,model,a.shape,flush=True);cnv.tl.infercnv(b,reference_key='group',reference_cat=refs,key_added=key,**spec['parameters']);C=b.obsm['X_'+key].tocsr();assert C.shape[0]==a.n_obs and np.isfinite(C.data).all();sparse.save_npz(D/(model+'_RNA_CNV.npz'),C);meta=b.uns[key];(D/(model+'_genomic_windows.json')).write_text(json.dumps(meta,indent=2,default=serial),encoding='utf-8');magnitude=np.asarray(abs(C).mean(1)).ravel();obs=b.obs.copy();obs['RNA_CNV_mean_absolute']=magnitude;obs.to_csv(D/(model+'_cell_descriptors.tsv'),sep='\t',encoding='utf-8');gmeans={}
  for group in group_order:
   mask=b.obs.group.eq(group).to_numpy()
   if not mask.any():continue
   gmeans[group]=np.asarray(C[mask].mean(0)).ravel();summaries.append({'dataset':item['dataset'],'sample':item['sample'],'reference_model':model,'group':group,'n_cells':int(mask.sum()),'median_mean_absolute_RNA_CNV':float(np.median(magnitude[mask])),'mean_mean_absolute_RNA_CNV':float(np.mean(magnitude[mask])),'high_complexity_cells':int(b.obs.loc[mask,'high_complexity_gt7000genes'].sum())})
  profiles[model]=gmeans
  for group in ['conservative_epithelial','immune_like','endothelial_like','fibroblast_like']:
   x=gmeans['strict_mesothelial'];y=gmeans[group];rho=float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else None;similarities.append({'dataset':item['dataset'],'sample':item['sample'],'reference_model':model,'comparison':'mesothelial_vs_'+group,'window_profile_correlation':rho,'biological_P':None})
  order=np.concatenate([np.flatnonzero(b.obs.group.eq(g)) for g in group_order]);assert len(order)==b.n_obs;arr=C[order].toarray();fig,ax=plt.subplots(figsize=(15,7));im=ax.imshow(arr,aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-.5,vmax=.5);ys=[];names=[];start=0
  for group in group_order:
   n=int(b.obs.group.eq(group).sum())
   if not n:continue
   ys.append(start+(n-1)/2);names.append(group+' (n='+str(n)+')');start+=n;ax.axhline(start-.5,color='#333333',lw=.4)
  ax.set_yticks(ys,names,fontsize=8);positions=sorted(meta['chr_pos'].items(),key=lambda x:x[1]);xt=[];xl=[]
  for i,(chrom,start) in enumerate(positions):
   end=positions[i+1][1] if i+1<len(positions) else C.shape[1];xt.append((start+end-1)/2);xl.append(chrom);ax.axvline(start-.5,color='#333333',lw=.3)
  ax.set_xticks(xt,xl,rotation=90,fontsize=8);ax.set_title(sample+' | '+model+' | RNA-derived chromosomal patterns',fontsize=12);ax.set_xlabel('Autosomal expression windows; not measured DNA copy numbers');fig.colorbar(im,ax=ax,label='Centered smoothed relative RNA expression',fraction=.025,pad=.02);fig.tight_layout();fig.savefig(D/(model+'_heatmap.png'),dpi=150);plt.close(fig);case['models'].append({'name':model,'reference':refs,'shape':list(C.shape),'finite_values':True,'plot_color_limit':[-.5,.5],'outside_color_range_fraction':float(np.count_nonzero(abs(C.data)>.5)/(C.shape[0]*C.shape[1])),'min':float(C.min()),'max':float(C.max())})
 (D/'CASE_COMPLETE.json').write_text(json.dumps(case,indent=2),encoding='utf-8');records.append(case);print('CNV_DONE',sample,flush=True)
pd.DataFrame(summaries).to_csv(root/'group_descriptors.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(similarities).to_csv(root/'mesothelial_profile_similarity.tsv',sep='\t',index=False,encoding='utf-8');out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'package_version':cnv.__version__,'cases':records,'all_pass':len(records)==len(spec['cases']),'interpretation':'auxiliary RNA-CNV identity evidence only; no categorical DNA diagnosis or biological P; review reference sensitivity and heatmaps before conclusions','visual_and_numeric_QA_pending':True};(O/'results/CNV_IDENTITY_COMPLETE_167.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
