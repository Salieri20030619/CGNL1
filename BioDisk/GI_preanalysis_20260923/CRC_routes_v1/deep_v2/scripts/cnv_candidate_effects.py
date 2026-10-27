import pathlib,json,numpy as np,pandas as pd
from scipy import sparse
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');B=O.parent;R=O/'results/cnv';X=sparse.load_npz(R/'counts_gene_by_cell.npz');sy=pd.read_csv(R/'gene_names.tsv',sep='\t').gene.values;m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);m.index=m.index.str.replace('-','.',regex=False);c=pd.read_csv(R/'cell_CNV_classification.tsv.gz',sep='\t',index_col=0);plan=json.loads((O/'PLAN.json').read_text());frozen=plan['frozen_genes'];rows=[];coverage=[];summ=[]
for quant in [95,99]:
 for exclude in [False,True]:
  mask=c['cnv_supported_'+str(quant)]&(~c.cluster.str.startswith('Tu10') if exclude else True);cc=c[mask];label='CNV'+str(quant)+('_noTu10' if exclude else '');eligible={}
  for (pt,site),z in cc.groupby(['patients','organs']):
   ix=m.index.get_indexer(z.index);assert (ix>=0).all();coverage.append({'subset':label,'patient':pt,'site':site,'cells':len(ix)})
   if len(ix)<20:continue
   pb=np.asarray(X[:,ix].sum(axis=1)).ravel();lib=pb.sum();lc=np.log2(1+pb/lib*1e6);cp=np.asarray(X[:,ix].multiply(1e4/m.iloc[ix].nCount_RNA.values).mean(axis=1)).ravel();eligible[(pt,site)]=(lc,cp,len(ix))
  pts=sorted({pt for pt,site in eligible if (pt,'CCT') in eligible and (pt,'LCT') in eligible});deltas=[]
  for pt in pts:
   pr,pm=eligible[(pt,'CCT')],eligible[(pt,'LCT')];d=pm[0]-pr[0];dm=pm[1]-pr[1];deltas.append(d)
   for i,g in enumerate(sy):
    if g in frozen:rows.append({'candidate':g,'subset':label,'patient':pt,'primary_log2CPM':pr[0][i],'liver_log2CPM':pm[0][i],'delta':d[i],'delta_mean_CP10K':dm[i],'primary_cells':pr[2],'liver_cells':pm[2]})
  ds=np.array(deltas);s=pd.DataFrame({'gene':sy,'subset':label,'n_pairs':len(pts),'median_delta':np.median(ds,axis=0),'positive':(ds>0).sum(axis=0),'negative':(ds<0).sum(axis=0)});summ.append(s)
r=pd.DataFrame(rows);r.to_csv(R/'candidate_patient_effects.tsv',sep='\t',index=False);s=pd.concat(summ);s.to_csv(R/'genomewide_CNV_subset_effects.tsv.gz',sep='\t',index=False);s[s.gene.isin(frozen)].to_csv(R/'candidate_CNV_summary.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(R/'subset_cell_coverage.tsv',sep='\t',index=False);(R/'EFFECTS_COMPLETE.json').write_text(json.dumps({'exit_code':0,'subsets':4,'candidates':len(frozen),'scope':'n4 paired effects; no inferential tests;95/99 reference thresholds and Tu10 exclusion are sensitivities'},indent=2));print(s[s.gene.isin(['CGNL1','CREB3L1','PCSK5','SLC9B2','BTNL8','TMEM131L'])].to_string(index=False))
