import pathlib,json,tarfile,gzip,numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'pair_validation';P=json.loads((O/'PLAN.json').read_text());lineage=json.loads((V/'00_docs/external_gating_plan.json').read_text())['lineage'];wanted=sorted(set(P['genes']+sum(P['programs'].values(),[])+P['gastric_proxy']+sum(lineage.values(),[])+['PTPRC','MUC1','CEACAM5','CEACAM6','MKI67','VIM','CDH1','CLDN18','TFF1','TFF2','TFF3']))
archive=O/'GSE308231_RAW.tar';rows=[];q=[];cells=[];features=[]
with tarfile.open(archive) as tar:
 names=tar.getnames();prefixes=sorted(n[:-len('_matrix.mtx.gz')] for n in names if n.endswith('_matrix.mtx.gz'));assert len(prefixes)==6
 for prefix in prefixes:
  def f(n):return gzip.GzipFile(fileobj=tar.extractfile(n))
  ft=pd.read_csv(f(prefix+'_features.tsv.gz'),sep='\t',header=None);sy=ft.iloc[:,1].astype(str).values;bc=pd.read_csv(f(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str).values
  with threadpool_limits(limits=2):X=io.mmread(f(prefix+'_matrix.mtx.gz')).tocsr()
  assert X.shape==(len(sy),len(bc));assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data));u=np.asarray(X.sum(axis=0)).ravel();ng=np.asarray((X>0).sum(axis=0)).ravel();mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(u,1)*100;keep=(u>=500)&(ng>=200)&(ng<=7000)&(mt<=25);X=X[:,keep];umi=u[keep];n=X.shape[1]
  group='PT' if '_Ca_' in prefix else 'PM';assert group=='PT' or '_F_' in prefix
  def detect(gs):return np.asarray((X[np.isin(sy,gs)]>0).sum(axis=0)).ravel()
  k={a:detect(gs) for a,gs in lineage.items()};my=k['myeloid_core']>=2;ly=(k['T']>=2)|(k['B']>=2);ep=k['epithelial']>=2;mes=k['mesothelial']>=2
  label=np.repeat('Other_unresolved',n).astype(object);label[ly&~my&~ep]='Lymphoid';label[my&~ly&~ep]='Myeloid_other';label[my&(k['macrophage_mature']>=2)&~ly&~ep]='Macrophage_enriched';label[ep&~my&~ly]='Epithelial_like';label[(k['fibroblast']>=2)&~ep&~my&~ly]='Fibroblast_like';label[(k['endothelial']>=2)&~ep&~my&~ly]='Endothelial_like';label[mes&~my&~ly]='Mesothelial_like';label[(my&ly)|(my&ep)|(ly&ep)]='Mixed_lineage'
  strict=(label=='Epithelial_like')&(detect(['EPCAM'])>0)&(k['epithelial']>=3)&(detect(['PTPRC'])==0)&(k['mesothelial']==0)
  gg=np.array(sorted(set(sy)&set(wanted)));Y=sparse.vstack([sparse.csr_matrix(X[np.flatnonzero(sy==g)].sum(axis=0)) for g in gg],format='csr');assert set(P['genes'])<=set(gg)
  features.append({'sample':prefix,'features':len(sy),'unique_symbols':len(set(sy)),'duplicate_symbols':';'.join(g for g in wanted if sum(sy==g)>1),'missing_requested':';'.join(sorted(set(wanted)-set(gg)))})
  cm=pd.DataFrame({'sample':prefix,'group':group,'barcode':bc[keep],'celltype':label,'strict_epithelial':strict,'umi':umi,'n_genes':ng[keep],'mt_percent':mt[keep]})
  for g in P['genes']:
   vals=np.asarray(Y[np.flatnonzero(gg==g)].sum(axis=0)).ravel();cm[g+'_counts']=vals;cm[g+'_CP10K']=vals/np.maximum(umi,1)*1e4
  cells.append(cm)
  groups={l:label==l for l in np.unique(label)};groups['Epithelial_strict']=strict;groups['All_QC']=np.ones(n,dtype=bool)
  for lab,mask in groups.items():
   ix=np.flatnonzero(mask)
   if len(ix)==0:continue
   cnt=np.asarray(Y[:,ix].sum(axis=1)).ravel();det=np.asarray((Y[:,ix]>0).sum(axis=1)).ravel();lib=umi[ix].sum();yn=Y[:,ix].multiply(1e4/np.maximum(umi[ix],1));mean_cp=np.asarray(yn.mean(axis=1)).ravel();yn=yn.tocsr();yn.data=np.log1p(yn.data);mean_ln=np.asarray(yn.mean(axis=1)).ravel()
   for j,g in enumerate(gg):rows.append({'sample':prefix,'group':group,'celltype':lab,'gene':g,'n_cells':len(ix),'counts':cnt[j],'library_size':lib,'n_detected':det[j],'fraction_detected':det[j]/len(ix),'log2CPM':np.log2(1+cnt[j]/max(lib,1)*1e6),'mean_CP10K':mean_cp[j],'mean_log1p_CP10K':mean_ln[j]})
  q.append({'sample':prefix,'group':group,'input_cells':len(bc),'qc_cells':n,'median_umi':float(np.median(umi)),'median_genes':float(np.median(ng[keep])),'epithelial_like':int((label=='Epithelial_like').sum()),'epithelial_strict':int(strict.sum()),'mesothelial_like':int((label=='Mesothelial_like').sum())});print(json.dumps(q[-1]),flush=True)
  # Keep full gene counts in epithelial cells for candidate-independent follow-up, if warranted.
  ix=np.flatnonzero(label=='Epithelial_like');sparse.save_npz(O/(prefix+'_epithelial_counts.npz'),X[:,ix]);pd.DataFrame({'gene_id':ft.iloc[:,0],'symbol':sy}).to_csv(O/(prefix+'_features.tsv.gz'),sep='\t',index=False);cm.iloc[ix].to_csv(O/(prefix+'_epithelial_metadata.tsv.gz'),sep='\t',index=False)
a=pd.DataFrame(rows);a.to_csv(O/'tissue_gene_compartment_values.tsv.gz',sep='\t',index=False);pd.DataFrame(q).to_csv(O/'tissue_sample_qc.tsv',sep='\t',index=False);pd.DataFrame(features).to_csv(O/'tissue_features_qc.tsv',sep='\t',index=False);pd.concat(cells).to_csv(O/'tissue_cell_metadata.tsv.gz',sep='\t',index=False)
summary=[]
for cutoff in [20,50]:
 for ct in ['Epithelial_like','Epithelial_strict','Mesothelial_like','Macrophage_enriched','All_QC']:
  for g in P['genes']:
   z=a[(a.celltype==ct)&(a.gene==g)&(a.n_cells>=cutoff)];pt=z[z.group=='PT'];pm=z[z.group=='PM'];r={'gene':g,'celltype':ct,'min_cells':cutoff,'PT_samples':len(pt),'PM_samples':len(pm)}
   for metric in ['log2CPM','fraction_detected','mean_CP10K','mean_log1p_CP10K']:
    r[metric+'_PT_median']=pt[metric].median();r[metric+'_PM_median']=pm[metric].median();r[metric+'_median_difference']=pm[metric].median()-pt[metric].median()
   if len(pt)>0 and len(pm)>0:
    dif=pm.log2CPM.values[:,None]-pt.log2CPM.values[None,:];r['cross_comparisons_lower']=int((dif<0).sum());r['cross_comparisons_total']=dif.size;r['all_PM_below_all_PT']=bool(pm.log2CPM.max()<pt.log2CPM.min())
   summary.append(r)
s=pd.DataFrame(summary);s.to_csv(O/'tissue_candidate_descriptive_comparison.tsv',sep='\t',index=False)
(O/'TISSUE_COMPLETE.json').write_text(json.dumps({'exit_code':0,'samples':len(q),'PT_samples':sum(r['group']=='PT' for r in q),'PM_samples':sum(r['group']=='PM' for r in q),'qc_cells':sum(r['qc_cells'] for r in q),'scope':'descriptive unpaired sample-level only; donor pairing not supplied; epithelial-like not confirmed malignant; no inferential P or FDR'},indent=2));print(s.to_string(index=False))
