import pathlib,json,numpy as np,pandas as pd
from scipy import sparse
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=O/'results/route1';R.mkdir(exist_ok=True);plan=json.loads((O/'metadata/ANALYSIS_PLAN.json').read_text());m=pd.read_csv(O/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);idx=m.index.str.replace('-','.',regex=False);assert idx.is_unique;m.index=idx
reader=pd.read_csv(O/'raw/nonimmune_counts.txt.gz',sep='\t',index_col=0,chunksize=128);groups=[];inds=[]
for (pt,site),z in m.groupby(['patients','organs']):
 assert pt in plan['clinical_map']
 for name,mask in [('Tumor_author',z.cluster.str.startswith('Tu')),('Tumor_no_Tu10',z.cluster.str.startswith('Tu')&~z.cluster.str.startswith('Tu10')),('Tumor_strict_QC',z.cluster.str.startswith('Tu')&~z.cluster.str.startswith('Tu10')&(z['percent.mt']<=25)&(z.nFeature_RNA>=200)),('Fibroblast',z.cluster.str.startswith('F')),('Endothelial',z.cluster.str.startswith('E'))]+[(str(c),z.cluster==c) for c in z.cluster.unique()]:
  ii=m.index.get_indexer(z.index[mask]);groups.append({'patient':pt,'site':site,'celltype':name,'n_cells':len(ii),'median_UMI':float(m.iloc[ii].nCount_RNA.median()) if len(ii) else np.nan,'median_MT':float(m.iloc[ii]['percent.mt'].median()) if len(ii) else np.nan});inds.append(ii)
g=pd.DataFrame(groups);g['group_id']=['G'+str(i) for i in range(len(g))];g.to_csv(R/'group_metadata.tsv',sep='\t',index=False);ri=np.concatenate(inds);ci=np.concatenate([np.repeat(i,len(v)) for i,v in enumerate(inds)]);B=sparse.csr_matrix((np.ones(len(ri)),(ri,ci)),shape=(len(m),len(g)));N=B.multiply(1e4/m.nCount_RNA.values[:,None]).tocsr();pb=[];mn=[];genes=[];celltotal=np.zeros(len(m));all_integer=True
for k,chunk in enumerate(reader):
 if k==0:assert set(chunk.columns)==set(m.index)
 x=chunk.loc[:,m.index].to_numpy(dtype=float);assert np.isfinite(x).all() and (x>=0).all();all_integer=all_integer and np.all(x==np.floor(x));celltotal+=x.sum(axis=0);pb.append((B.T@x.T).T);mn.append((N.T@x.T).T/np.maximum(g.n_cells.values,1));genes.extend(chunk.index.astype(str));
 if k%40==0:print('READ_GENES',len(genes),flush=True)
assert all_integer;check=np.max(np.abs(celltotal-m.nCount_RNA.values));print('MAX_COUNT_DIFF',check,flush=True);assert check<1e-6
P=np.vstack(pb);M=np.vstack(mn);sy=np.asarray(genes);assert len(set(genes))==len(genes);lib=P.sum(axis=0);cpm=P/np.maximum(lib,1)*1e6;L=np.log2(1+cpm);np.savez_compressed(R/'aggregates.npz',counts=P,mean_CP10K=M,log2CPM=L,genes=sy,group_ids=g.group_id.values.astype(str));pd.DataFrame(L,index=sy,columns=g.group_id).to_csv(R/'all_gene_log2CPM.tsv.gz',sep='\t',index_label='gene')
patientrows=[];summary=[];programrows=[]
for celltype in g.celltype.unique():
 z=g[(g.celltype==celltype)&(g.n_cells>=20)];pts=sorted(set(z[z.site=='CCT'].patient)&set(z[z.site=='LCT'].patient))
 if not pts:continue
 ds=[];dm=[];det=[]
 for pt in pts:
  i=z[(z.patient==pt)&(z.site=='CCT')].index[0];j=z[(z.patient==pt)&(z.site=='LCT')].index[0];d=L[:,j]-L[:,i];u=M[:,j]-M[:,i];ds.append(d);dm.append(u);det.append((cpm[:,i]>=1)|(cpm[:,j]>=1))
  for name,gset in plan['programs'].items():
   mask=np.isin(sy,gset);assert mask.sum()>=.8*len(gset);programrows.append({'patient':pt,'celltype':celltype,'program':name,'primary_score':L[mask,i].mean(),'liver_score':L[mask,j].mean(),'delta':d[mask].mean(),'n_primary':g.loc[i,'n_cells'],'n_liver':g.loc[j,'n_cells']})
  if celltype in ['Tumor_author','Tumor_no_Tu10','Tumor_strict_QC','Fibroblast']:
   patientrows.append(pd.DataFrame({'gene':sy,'patient':pt,'celltype':celltype,'delta_log2CPM':d,'delta_mean_CP10K':u}))
 ds=np.array(ds);dm=np.array(dm);det=np.array(det);ss=pd.DataFrame({'gene':sy,'celltype':celltype,'n_patients':len(pts),'median_delta_log2CPM':np.median(ds,axis=0),'mean_delta_log2CPM':ds.mean(axis=0),'positive':(ds>0).sum(axis=0),'negative':(ds<0).sum(axis=0),'expressed_patient_pairs':det.sum(axis=0),'median_delta_mean_CP10K':np.median(dm,axis=0),'min_delta':ds.min(axis=0),'max_delta':ds.max(axis=0)});summary.append(ss)
s=pd.concat(summary);s.to_csv(R/'genomewide_descriptive_effects.tsv.gz',sep='\t',index=False);pd.concat(patientrows).to_csv(R/'gene_patient_effects.tsv.gz',sep='\t',index=False);pp=pd.DataFrame(programrows);pp.to_csv(R/'program_patient_effects.tsv',sep='\t',index=False);ps=pp.groupby(['celltype','program']).agg(n=('patient','nunique'),median_delta=('delta','median'),mean_delta=('delta','mean'),positive=('delta',lambda x:int((x>0).sum())),negative=('delta',lambda x:int((x<0).sum()))).reset_index();ps.to_csv(R/'program_summary.tsv',sep='\t',index=False)
sel=s[(s.celltype=='Tumor_author')&(s.n_patients>=3)&(s.expressed_patient_pairs>=3)&(s.median_delta_log2CPM.abs()>=1)&(s[['positive','negative']].max(axis=1)>=np.ceil(.8*s.n_patients))&~s.gene.str.startswith(('MT-','RPL','RPS'))].copy();sel['abs_effect']=sel.median_delta_log2CPM.abs();sel.sort_values(['n_patients','abs_effect'],ascending=False).to_csv(R/'descriptive_screen.tsv',sep='\t',index=False)
coverage=[]
for ct,z in g.groupby('celltype'):
 zz=z[z.n_cells>=20];pts=sorted(set(zz[zz.site=='CCT'].patient)&set(zz[zz.site=='LCT'].patient));coverage.append({'celltype':ct,'paired_patients':len(pts),'patients':';'.join(pts)})
pd.DataFrame(coverage).to_csv(R/'paired_coverage.tsv',sep='\t',index=False);(R/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'genes':len(sy),'cells':len(m),'patients':m.patients.nunique(),'integer_counts':bool(all_integer),'max_cell_count_mismatch':float(check),'scope':'Author tumor labels; no CNV confirmation. All 5 donors have colon primaries and preoperative chemotherapy; descriptive only, no P/FDR.'},indent=2));print(pd.DataFrame(coverage).to_string(index=False));print(ps[ps.celltype.isin(['Tumor_author','Tumor_strict_QC'])].to_string(index=False));print(sel.sort_values('abs_effect',ascending=False).head(15).to_string(index=False))
