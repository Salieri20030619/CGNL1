from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib,datetime,gzip,tarfile,io,itertools
import numpy as np,pandas as pd
from scipy import stats,io as sio,sparse
import statsmodels.api as sm
from threadpoolctl import threadpool_limits
O=Path.cwd();B=O.parent;V=B/'GC_validation_v2';spec=json.loads((O/'ANALYSIS_SPECIFICATION.json').read_text());assert spec['target']=='CGNL1'
panels={'epithelial':['EPCAM','KRT8','KRT18','KRT19'],'stroma':['COL1A1','COL1A2','DCN','LUM'],'myeloid':['LYZ','LST1','TYROBP','FCER1G'],'adipocyte':['ADIPOQ','PLIN1','CIDEC','CIDEA'],'endothelial':['PECAM1','VWF','KDR','EMCN'],'mesothelial':['MSLN','CALB2','WT1','UPK3B']}
implementation={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_target_effects_not_yet_viewed':True,'bulk_panels':panels,'scores':'mean of original TMMlogCPM for fixed genes; >=3of4 required','singlecell_gates':'reuse frozentarget-excludedGSE228598QC/labels;GSE234129authorgroupsmappedconservatively','primary_test':'pairedWilcoxon exact when no zeros/ties; alsoenumerate conditionalrankdistribution if needed','M2':'source proxies are nuisance sensitivity, not causal mediation'}
(O/'IMPLEMENTATION_FREEZE.json').write_text(json.dumps(implementation,indent=2),encoding='utf-8')
inputs=[B/'03_processed/GC314812_TMM_logCPM.tsv.gz',B/'02_raw/GSE314812_counts.tsv.gz',B/'01_metadata/GSE314812_paired_design.tsv',V/'02_raw/GSE228598_RAW.tar',B/'02_raw/GSE234129_count_matrix.mtx.gz'];manifest=[]
for p in inputs:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for buf in iter(lambda:f.read(8*1024*1024),b''):h.update(buf)
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
old=json.loads((V/'logs/GSE228598_complete.json').read_text());assert next(x['sha256'] for x in manifest if 'GSE228598_RAW' in x['path'])==old['sha256']
(O/'metadata/REUSED_INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
dm=pd.read_csv(B/'01_metadata/GSE314812_paired_design.tsv',sep='\t');x=pd.read_csv(B/'03_processed/GC314812_TMM_logCPM.tsv.gz',sep='\t');q=x[x.symbol=='CGNL1'];assert len(q)==1
meta=dm.copy();meta['CGNL1']=q[meta['sample']].to_numpy().ravel();coverage=[]
for name,genes in panels.items():
 z=x[x.symbol.isin(genes)].groupby('symbol')[meta['sample']].mean();coverage.append({'panel':name,'genes_measured':len(z),'genes':'|'.join(z.index),'expected':'|'.join(genes)});meta[name]=z.mean(0).to_numpy() if len(z)>=3 else np.nan
pd.DataFrame(coverage).to_csv(O/'metadata/bulk_panel_coverage.tsv',sep='\t',index=False);meta.to_csv(O/'processed/bulk_sample_measurements.tsv',sep='\t',index=False)
delta=[]
for patient,g in meta.groupby('patient'):
 assert len(g)==2;pm=g[g.site=='GCPM'].iloc[0];pt=g[g.site=='GCPT'].iloc[0];row={'patient':patient,'PM_sample':pm['sample'],'PT_sample':pt['sample'],'PM_preservation':pm['preservation method'],'PT_preservation':pt['preservation method'],'preservation_delta':int(pm['preservation method']=='FFPE')-int(pt['preservation method']=='FFPE'),'primary':float(pt.CGNL1),'metastasis':float(pm.CGNL1)}
 for gene in ['CGNL1']+list(panels):row[gene]=float(pm[gene]-pt[gene])
 delta.append(row)
d=pd.DataFrame(delta);assert len(d)==17;d.to_csv(O/'results/bulk_patient_deltas.tsv',sep='\t',index=False)
values=d.CGNL1.to_numpy();nz=values[values!=0];ranks=stats.rankdata(abs(nz));dist=np.array([0.])
for rank in ranks:dist=np.r_[dist,dist+rank]
obs=min(ranks[nz>0].sum(),ranks[nz<0].sum());pw=min(1,2*np.mean(dist<=obs+1e-10))
result={'n':len(d),'mean_delta':float(values.mean()),'median_delta':float(np.median(values)),'positive':int((values>0).sum()),'negative':int((values<0).sum()),'P_exact_signed_rank':float(pw),'FDR':'not_required_single_prespecified_target','LOO_min':float(min(np.delete(values,i).mean() for i in range(len(values)))),'LOO_max':float(max(np.delete(values,i).mean() for i in range(len(values)))),'same_preservation_n':int((d.preservation_delta==0).sum()),'same_preservation_mean':float(d.loc[d.preservation_delta==0,'CGNL1'].mean()),'same_preservation_positive':int((d.loc[d.preservation_delta==0,'CGNL1']>0).sum())}
models=[]
for name,covs in [('M0',[]),('M1',['preservation_delta']),('M2',['preservation_delta']+list(panels))]:
 Z=np.column_stack([np.ones(len(d))]+[d[c].to_numpy() for c in covs]);rank=np.linalg.matrix_rank(Z);df=len(d)-rank
 if not np.isfinite(Z).all() or rank<Z.shape[1] or df<5:models.append({'model':name,'status':'not_estimable','n':len(d),'rank':rank,'df':df,'P':np.nan});continue
 fit=sm.OLS(values,Z).fit(cov_type='HC3',use_t=True);ci=fit.conf_int()[0];models.append({'model':name,'status':'estimated','n':len(d),'df':df,'effect':float(fit.params[0]),'P':float(fit.pvalues[0]),'CI_low':float(ci[0]),'CI_high':float(ci[1]),'condition':float(np.linalg.cond(Z))})
mod=pd.DataFrame(models);valid=mod.P.notna();mod['FDR_three_models']=np.nan;mod.loc[valid,'FDR_three_models']=stats.false_discovery_control(np.r_[mod.loc[valid,'P'],np.ones(3-valid.sum())])[:valid.sum()];mod.to_csv(O/'results/bulk_models.tsv',sep='\t',index=False)
# Independent normalization sensitivity using full-library expected-count denominator.
raw=pd.read_csv(B/'02_raw/GSE314812_counts.tsv.gz',sep='\t');hit=raw[raw.gene_id.astype(str).str.split('.').str[0]=='ENSG00000128849'];assert len(hit)==1;samples=dm['sample'].tolist();den=raw[samples].sum();cp=np.log2(hit[samples].iloc[0]/den*1e6+.5);rr=[]
for _,r in d.iterrows():rr.append({'patient':r.patient,'naive_logCPM_delta':float(cp[r.PM_sample]-cp[r.PT_sample])})
pd.DataFrame(rr).to_csv(O/'results/bulk_normalization_sensitivity.tsv',sep='\t',index=False);result['naive_CPM_mean_delta']=float(np.mean([r['naive_logCPM_delta'] for r in rr]));result['naive_CPM_positive']=sum(r['naive_logCPM_delta']>0 for r in rr)
(O/'results/BULK_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('BULK',json.dumps(result),mod.to_string(index=False),flush=True)
# GSE234129: actual matrix and author annotation, not inferred malignant labels.
cm=pd.read_csv(B/'GC_validation_v2/03_processed/GSE234129_annotation_qc.tsv.gz',sep='\t');bc=pd.read_csv(B/'02_raw/GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str);features=pd.read_csv(B/'02_raw/GSE234129_features.tsv.gz',sep='\t',header=None);genes=features.iloc[:,1 if features.shape[1]>1 else 0].astype(str).to_numpy()
with threadpool_limits(limits=4):X=sio.mmread(B/'02_raw/GSE234129_count_matrix.mtx.gz').tocsr()
ix=pd.Index(bc).get_indexer(cm.cell_barcodes.astype(str));assert (ix>=0).all();X=X[:,ix];assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data));tg=np.flatnonzero(genes=='CGNL1');assert len(tg)>0;target=np.asarray(X[tg].sum(0)).ravel();umi=np.asarray(X.sum(0)).ravel();cellrows=[]
for (patient,sample,celltype),ids in cm.groupby(['patient','sample','celltype']).groups.items():
 ids=np.array(list(ids));cellrows.append({'dataset':'GSE234129','patient':patient,'sample':sample,'celltype':celltype,'n_cells':len(ids),'target_counts':float(target[ids].sum()),'total_UMI':float(umi[ids].sum()),'detected':int((target[ids]>0).sum()),'fraction_detected':float((target[ids]>0).mean()),'ln1pCP10k':float(np.log1p(target[ids].sum()/max(umi[ids].sum(),1)*1e4))})
pd.DataFrame(cellrows).to_csv(O/'results/GSE234129_author_source.tsv',sep='\t',index=False);print('234_SOURCE',pd.DataFrame(cellrows).groupby('celltype')[['n_cells','target_counts']].sum().to_string(),flush=True);del X
# GSE228598: read raw sample matrices, join previously frozen target-blind labels by barcode.
rows=[];qc=[]
with tarfile.open(V/'02_raw/GSE228598_RAW.tar') as tf:
 prefixes=sorted(n[:-len('_matrix.mtx.gz')] for n in tf.getnames() if n.endswith('_matrix.mtx.gz'))
 for prefix in prefixes:
  def read(n):return gzip.GzipFile(fileobj=tf.extractfile(n))
  ft=pd.read_csv(read(prefix+'_features.tsv.gz'),sep='\t',header=None);genes=ft[1].astype(str).to_numpy();bc=pd.read_csv(read(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str)
  with threadpool_limits(limits=2):X=sio.mmread(read(prefix+'_matrix.mtx.gz')).tocsr()
  m=pd.read_csv(V/'03_processed'/(prefix+'_gating_refined.tsv.gz'),sep='\t');ix=pd.Index(bc).get_indexer(m.barcode.astype(str));assert (ix>=0).all();X=X[:,ix];tg=np.flatnonzero(genes=='CGNL1');assert len(tg)>0;t=np.asarray(X[tg].sum(0)).ravel();u=np.asarray(X.sum(0)).ravel();assert np.array_equal(u,m.umi)
  for lab,ids in m.groupby('label').groups.items():
   ids=np.array(list(ids));rows.append({'dataset':'GSE228598','patient':m.patient.iloc[0],'sample':prefix,'sample_type':m.sample_type.iloc[0],'celltype':lab,'n_cells':len(ids),'target_counts':float(t[ids].sum()),'total_UMI':float(u[ids].sum()),'detected':int((t[ids]>0).sum()),'fraction_detected':float((t[ids]>0).mean()),'ln1pCP10k':float(np.log1p(t[ids].sum()/max(u[ids].sum(),1)*1e4))})
  qc.append({'sample':prefix,'QC_cells':len(m),'target_total':float(t.sum())});print('228_SOURCE',prefix,'target_UMI',int(t.sum()),flush=True)
pd.DataFrame(rows).to_csv(O/'results/GSE228598_source.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(O/'metadata/GSE228598_reuse_QA.tsv',sep='\t',index=False)
summ=[]
for dataset,rs in [('GSE234129',cellrows),('GSE228598',rows)]:
 rr=pd.DataFrame(rs)
 for threshold in [20,50]:
  eligible=rr[rr.n_cells>=threshold];pat=eligible.groupby(['patient','celltype'])[['fraction_detected','ln1pCP10k']].mean().reset_index()
  for ct,g in pat.groupby('celltype'):summ.append({'dataset':dataset,'min_cells':threshold,'celltype':ct,'patients':g.patient.nunique(),'median_detection':float(g.fraction_detected.median()),'median_expression':float(g.ln1pCP10k.median())})
pd.DataFrame(summ).to_csv(O/'results/source_patient_summary.tsv',sep='\t',index=False)
(O/'BULK_SOURCE_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bulk_patients':17,'independent_source_patients':len(qc),'raw_source_counts_reconstructed':True,'no_cell_level_P':True},indent=2),encoding='utf-8');print('SOURCE_SUMMARY',pd.DataFrame(summ).to_string(index=False))
