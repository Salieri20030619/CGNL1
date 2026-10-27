from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,gzip,tarfile,io
import numpy as np,pandas as pd
from scipy import io as sio,stats
import statsmodels.api as sm
from threadpoolctl import threadpool_limits
O=Path.cwd();B=O.parent;V=B/'GC_validation_v2';S=json.loads((O/'MESOTHELIAL_IDENTITY_SPEC.json').read_text());ip=S['identity_panels'];sp=S['state_panels'];needed=['CGNL1']+sum(ip.values(),[])+sum(sp.values(),[]);rows=[];coverage=[]
with tarfile.open(V/'02_raw/GSE228598_RAW.tar') as tf:
 prefixes=sorted(n[:-len('_matrix.mtx.gz')] for n in tf.getnames() if n.endswith('_matrix.mtx.gz'))
 for prefix in prefixes:
  def read(n):return gzip.GzipFile(fileobj=tf.extractfile(n))
  ft=pd.read_csv(read(prefix+'_features.tsv.gz'),sep='\t',header=None);genes=ft[1].astype(str).to_numpy();bc=pd.read_csv(read(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str)
  with threadpool_limits(limits=2):X=sio.mmread(read(prefix+'_matrix.mtx.gz')).tocsr()
  m=pd.read_csv(V/'03_processed'/(prefix+'_gating_refined.tsv.gz'),sep='\t');ix=pd.Index(bc).get_indexer(m.barcode.astype(str));assert (ix>=0).all();X=X[:,ix];u=np.asarray(X.sum(0)).ravel();assert np.array_equal(u,m.umi);counts={g:np.asarray(X[np.flatnonzero(genes==g)].sum(0)).ravel() for g in set(needed) if g in genes}
  def det(gs):return sum((counts[g]>0).astype(int) for g in gs if g in counts)
  my=det(ip['myeloid'])>=2;ly=m.label.eq('Lymphoid').to_numpy();mes=(det(ip['mesothelial_specific'])>=1)&(det(ip['mesothelial_support'])>=1)&~my&~ly;epi=(counts['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)&(det(ip['mesothelial_specific'])==0)&~my&~ly;cancer=det(ip['cancer_support'])>=2
  masks={'strict_mesothelial':mes,'conservative_epithelial':epi,'legacy_mesothelial':m.label.eq('Mesothelial_like').to_numpy(),'macrophage':m.label.eq('Macrophage_enriched').to_numpy(),'strict_mesothelial_low_cancer_proxy':mes&~cancer}
  cells=m[['barcode','patient','sample']].copy();cells['strict_mesothelial']=mes;cells['conservative_epithelial']=epi;cells['cancer_proxy_ge2']=cancer;cells['CGNL1_count']=counts['CGNL1'];cells.to_csv(O/'processed'/(prefix+'_identity_validation.tsv.gz'),sep='\t',index=False)
  for group,mask in masks.items():
   n=int(mask.sum());target=counts['CGNL1'][mask];total=float(u[mask].sum());row={'patient':m.patient.iloc[0],'sample':prefix,'group':group,'n_cells':n,'total_UMI':total,'mean_logUMI':float(np.log1p(u[mask]).mean()) if n else np.nan,'CGNL1_count':float(target.sum()),'CGNL1_detected':int((target>0).sum()),'CGNL1_detection':float((target>0).mean()) if n else np.nan,'CGNL1':float(np.log1p(target.sum()/max(total,1)*1e4)) if n else np.nan,'cancer_proxy_cells':int((mask&cancer).sum())}
   for name,gs in dict(sp,mesothelial_identity=ip['mesothelial_specific']+ip['mesothelial_support']).items():row[name]=float(np.mean([np.log1p(counts[g][mask].sum()/max(total,1)*1e4) for g in gs if g in counts])) if n else np.nan
   rows.append(row)
  coverage.append({'patient':m.patient.iloc[0],'strict_mesothelial':int(mes.sum()),'conservative_epithelial':int(epi.sum()),'legacy_mesothelial':int(masks['legacy_mesothelial'].sum()),'strict_low_cancer_proxy':int(masks['strict_mesothelial_low_cancer_proxy'].sum()),'strict_fraction_also_legacy':float(masks['legacy_mesothelial'][mes].mean()) if mes.any() else np.nan});print('IDENTITY',prefix,int(mes.sum()),int(epi.sum()),flush=True)
q=pd.DataFrame(rows);q.to_csv(O/'results/identity_source_profiles.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'results/identity_coverage.tsv',sep='\t',index=False);summary=[];contr=[];states=[]
for cutoff in [20,50]:
 z=q[q.n_cells>=cutoff]
 for group,g in z.groupby('group'):summary.append({'min_cells':cutoff,'group':group,'patients':g.patient.nunique(),'median_CGNL1':float(g.CGNL1.median()),'median_detection':float(g.CGNL1_detection.median()),'patients_detected':int((g.CGNL1_count>0).sum())})
 for other in ['conservative_epithelial','macrophage']:
  p=z[z.group.isin(['strict_mesothelial',other])].pivot(index='patient',columns='group',values='CGNL1').dropna();d=p.strict_mesothelial-p[other] if len(p) and other in p else pd.Series(dtype=float);P=float(stats.wilcoxon(d,method='auto').pvalue) if len(d)>=10 and (d!=0).any() else np.nan;contr.append({'min_cells':cutoff,'contrast':'strict_mesothelial-minus-'+other,'n':len(d),'mean':float(d.mean()) if len(d) else np.nan,'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P':P})
 mes=z[z.group=='strict_mesothelial']
 for panel in sp:
  if len(mes)<10:states.append({'min_cells':cutoff,'panel':panel,'model':'both','n':len(mes),'status':'descriptive_only_small_n','P':np.nan});continue
  r=stats.spearmanr(mes.CGNL1,mes[panel]);states.append({'min_cells':cutoff,'panel':panel,'model':'Spearman','n':len(mes),'status':'estimated','effect':float(r.statistic),'P':float(r.pvalue)})
  X=np.column_stack([np.ones(len(mes)),mes.CGNL1,mes.mean_logUMI,mes.mesothelial_identity]);rank=np.linalg.matrix_rank(X)
  if rank==X.shape[1] and len(mes)-rank>=5:
   fit=sm.OLS(mes[panel],X).fit(cov_type='HC3',use_t=True);states.append({'min_cells':cutoff,'panel':panel,'model':'HC3_depth_identity','n':len(mes),'status':'estimated','effect':float(fit.params[1]),'P':float(fit.pvalues[1])})
c=pd.DataFrame(contr);c['FDR']=np.nan
for cutoff,ix in c.groupby('min_cells').groups.items():
 good=[i for i in ix if np.isfinite(c.loc[i,'P'])]
 if good:c.loc[good,'FDR']=stats.false_discovery_control(np.r_[c.loc[good,'P'],np.ones(2-len(good))])[:len(good)]
s=pd.DataFrame(states);s['FDR']=np.nan
for key,ix in s.groupby(['min_cells','model']).groups.items():
 good=[i for i in ix if np.isfinite(s.loc[i,'P'])]
 if good:s.loc[good,'FDR']=stats.false_discovery_control(np.r_[s.loc[good,'P'],np.ones(3-len(good))])[:len(good)]
pd.DataFrame(summary).to_csv(O/'results/identity_source_summary.tsv',sep='\t',index=False);c.to_csv(O/'results/identity_source_contrasts.tsv',sep='\t',index=False);s.to_csv(O/'results/mesothelial_state_models.tsv',sep='\t',index=False);print('SUMMARY',pd.DataFrame(summary).to_string(index=False));print('CONTRASTS',c.to_string(index=False));print('STATES',s.to_string(index=False));(O/'IDENTITY_VALIDATION_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'patients':q.patient.nunique(),'new_gates_target_excluded':True,'no_PC_status_inferred_from_sample_type':True},indent=2),encoding='utf-8')
