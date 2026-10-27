from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,gzip,datetime,hashlib
import numpy as np,pandas as pd
from scipy import sparse,io,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'reassessment_20261008_v1';O=R/'results';spec=json.loads((R/'REASSESSMENT_SPEC.json').read_text());old=B/'review_consensus_retest_20261002_v1'
programs=spec['R2_lineage_state']['programs'];hyp=json.loads((B/'v3_pilot_20261006/plan_v3.json').read_text())['panels']['hypoxia'];needed={'CGNL1'}|set(hyp)|set(sum(programs.values(),[]));rows=[]
reg=json.loads((old/'CACHE_INDEX.json').read_text());anns={d:pd.read_csv(old/'results'/f'{d}_accepted_annotation.tsv.gz',sep='\t') for d in sorted({r['dataset'] for r in reg})}
for rec in reg:
 print('LINEAGE_START',rec['sample'],flush=True);stem=rec['source_matrix'].replace('matrix.mtx.gz','');genes=pd.read_csv(stem+'features.tsv.gz',sep='\t',header=None)[1].astype(str).to_numpy();bc=pd.read_csv(stem+'barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
 with gzip.open(rec['source_matrix'],'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 a=anns[rec['dataset']];a=a[a['sample'].eq(rec['sample'])].copy();ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all();X=X[:,ix];umi=np.asarray(X.sum(0)).ravel();assert np.array_equal(umi,a.UMI_QC)
 for label in spec['R2_lineage_state']['lineages']:
  mask=a.accepted_major_lineage.eq(label).to_numpy();n=int(mask.sum())
  row=dict(dataset=rec['dataset'],sample=rec['sample'],patient=rec['donor'],lineage=label,cells=n,eligible=n>=20)
  if n:
   total=np.asarray(X[:,mask].sum(1)).ravel();den=total.sum();vals={g:float(np.log1p(total[genes==g].sum()/max(den,1)*1e4)) for g in needed if np.any(genes==g)}
   row.update(total_UMI=float(den),mean_logUMI=float(np.log1p(umi[mask]).mean()),CGNL1=vals.get('CGNL1',np.nan),CGNL1_count=float(total[genes=='CGNL1'].sum()),hypoxia=float(np.mean([vals[g] for g in hyp if g in vals])))
   for k,gs in programs.items():
    available=[vals[g] for g in gs if g in vals];row[k]=float(np.mean(available)) if len(available)>=.8*len(gs) else np.nan;row[k+'_covered']=len(available)
  rows.append(row)
 pd.DataFrame(rows).to_csv(O/'R2_library_lineage_profiles.tsv',sep='\t',index=False);print('LINEAGE_DONE',rec['sample'],flush=True)
lib=pd.DataFrame(rows);pt=lib[lib.eligible].groupby(['dataset','patient','lineage'],as_index=False).mean(numeric_only=True);pt.to_csv(O/'R2_patient_lineage_profiles.tsv',sep='\t',index=False)
models=[];pred=[];status=[]
for scope in ['pooled','GSE245552','GSE285991']:
 for label in spec['R2_lineage_state']['lineages']:
  v=pt[pt.lineage.eq(label)];v=v if scope=='pooled' else v[v.dataset.eq(scope)];v=v.reset_index(drop=True)
  for program in programs:
   vv=v.dropna(subset=['CGNL1','mean_logUMI','hypoxia',program]).reset_index(drop=True);n=len(vv)
   if n<4:status.append(dict(scope=scope,lineage=label,program=program,n=n,status='insufficient_patients'));continue
   base=pd.DataFrame({'const':np.ones(n),'depth':vv.mean_logUMI,'hypoxia':vv.hypoxia});base=pd.concat([base,pd.get_dummies(vv.dataset,drop_first=True,dtype=float)],axis=1);F=base.copy();F['CGNL1']=vv.CGNL1;y=vv[program]
   if np.linalg.matrix_rank(F)<F.shape[1] or n-F.shape[1]<3 or y.std()<1e-12:status.append(dict(scope=scope,lineage=label,program=program,n=n,status='rank_or_residual_df_or_no_variation'));continue
   status.append(dict(scope=scope,lineage=label,program=program,n=n,status='estimated_exploratory'))
   for kind in ['OLS','HC3']:
    m=sm.OLS(y,F).fit(cov_type='nonrobust' if kind=='OLS' else 'HC3',use_t=True);ci=m.conf_int().loc['CGNL1'];models.append(dict(scope=scope,lineage=label,program=program,variance=kind,n=n,beta=float(m.params.CGNL1)*vv.CGNL1.std()/y.std(),slope=float(m.params.CGNL1),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(m.pvalues.CGNL1),family='R2_'+scope+'_'+program+'_'+kind,evidence='exploratory_transcriptional_state_not_protein_function'))
   if n-F.shape[1]>=4:
    for i in range(n):
     keep=np.arange(n)!=i;f0=np.linalg.lstsq(base.to_numpy()[keep],y.to_numpy()[keep],rcond=None)[0];f1=np.linalg.lstsq(F.to_numpy()[keep],y.to_numpy()[keep],rcond=None)[0]
     pred.append(dict(scope=scope,lineage=label,program=program,patient=vv.patient.iloc[i],observed=float(y.iloc[i]),baseline=float(base.to_numpy()[i]@f0),with_CGNL1=float(F.to_numpy()[i]@f1)))
t=pd.DataFrame(models);t['FDR']=np.nan
for f,ix in t.groupby('family').groups.items():t.loc[ix,'FDR']=multipletests(t.loc[ix,'P'],method='fdr_bh')[1]
t.to_csv(O/'R2_state_models.tsv',sep='\t',index=False);pd.DataFrame(status).to_csv(O/'R2_model_eligibility.tsv',sep='\t',index=False);pr=pd.DataFrame(pred);pr.to_csv(O/'R2_heldout_predictions.tsv',sep='\t',index=False)
perf=[]
for key,g in pr.groupby(['scope','lineage','program']):
 b=(g.observed-g.baseline)**2;a=(g.observed-g.with_CGNL1)**2;perf.append(dict(zip(['scope','lineage','program'],key),n=len(g),MSE_baseline=float(b.mean()),MSE_with_CGNL1=float(a.mean()),relative_MSE_improvement=float(1-a.mean()/b.mean()),patients_improved=int((a<b).sum())))
pd.DataFrame(perf).to_csv(O/'R2_incremental_performance.tsv',sep='\t',index=False)

# R3: disease-site contrasts within fixed cell sources and matched patients.
agg=np.load(B/'results/route1/aggregates.npz',allow_pickle=True);gs=agg['genes'].astype(str);gi=np.flatnonzero(gs=='CGNL1');assert len(gi)==1;metadata=pd.read_csv(B/'results/route1/group_metadata.tsv',sep='\t');assert np.array_equal(metadata.group_id,agg['group_ids'])
counts=agg['counts'];L=np.log2(1+counts[gi[0]]/np.maximum(counts.sum(0),1)*1e6);pairs=[]
for celltype in spec['R3_disease_contrast']['fixed_groups']:
 v=metadata[metadata.celltype.eq(celltype)&metadata.n_cells.ge(20)]
 for patient in sorted(set(v.loc[v.site.eq('CCT'),'patient'])&set(v.loc[v.site.eq('LCT'),'patient'])):
  i=v[v.patient.eq(patient)&v.site.eq('CCT')].index[0];j=v[v.patient.eq(patient)&v.site.eq('LCT')].index[0]
  pairs.append(dict(group=celltype,patient=patient,primary=float(L[i]),liver=float(L[j]),delta=float(L[j]-L[i]),primary_cells=int(metadata.loc[i,'n_cells']),liver_cells=int(metadata.loc[j,'n_cells']),family='R3_main_three_sources'))
# Recompute frozen CNV95/99 exclusion sensitivities from raw integer count cache.
cdir=B/'deep_v2/results/cnv';X=sparse.load_npz(cdir/'counts_gene_by_cell.npz');gn=pd.read_csv(cdir/'gene_names.tsv',sep='\t').gene.to_numpy();ixgene=np.flatnonzero(gn=='CGNL1');assert len(ixgene)==1
m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);m.index=m.index.str.replace('-','.',regex=False);cnv=pd.read_csv(cdir/'cell_CNV_classification.tsv.gz',sep='\t',index_col=0);assert X.shape[1]==len(m)
for q in [95,99]:
 subset=cnv[cnv['cnv_supported_'+str(q)]&~cnv.cluster.str.startswith('Tu10')];eligible={}
 for (patient,site),a in subset.groupby(['patients','organs']):
  if len(a)<20:continue
  ix=m.index.get_indexer(a.index);assert (ix>=0).all();den=float(X[:,ix].sum());val=float(X[ixgene[0],ix].sum());eligible[(patient,site)]=(np.log2(1+val/den*1e6),len(ix))
 for patient in sorted({p for p,s in eligible if (p,'CCT') in eligible and (p,'LCT') in eligible}):
  a,n1=eligible[(patient,'CCT')];b,n2=eligible[(patient,'LCT')];pairs.append(dict(group='CNV'+str(q)+'_noTu10',patient=patient,primary=a,liver=b,delta=b-a,primary_cells=n1,liver_cells=n2,family='R3_CNV_two_sensitivities'))
pa=pd.DataFrame(pairs);pa.to_csv(O/'R3_paired_disease_effects.tsv',sep='\t',index=False);summary=[]
for group,g in pa.groupby('group'):
 v=g.delta.to_numpy();nz=v[np.abs(v)>1e-12];pv=stats.binomtest(int((nz>0).sum()),len(nz),.5).pvalue if len(nz) else np.nan
 summary.append(dict(group=group,n=len(g),positive=int((v>0).sum()),negative=int((v<0).sum()),mean_delta=float(v.mean()),median_delta=float(np.median(v)),P=float(pv),minimum_two_sided_sign_P=float(2/2**len(nz)) if len(nz) else np.nan,family=g.family.iloc[0],evidence='small_sample_exploratory_paired_site_contrast'))
t=pd.DataFrame(summary);t['FDR']=np.nan
for f,ix in t.groupby('family').groups.items():t.loc[ix,'FDR']=multipletests(t.loc[ix,'P'],method='fdr_bh')[1]
t.to_csv(O/'R3_paired_summary.tsv',sep='\t',index=False)
bulks=[]
for cohort,rel in [('bulk_validation','deep_v2/results/bulk_validation'),('second_validation','deep_v2/results/second_validation')]:
 t=pd.read_csv(B/rel/'candidate_validation_all_models.tsv',sep='\t');t=t[t.candidate.eq('CGNL1')].copy();t['source']=str(B/rel);t['cohort']=cohort;bulks.append(t)
pd.concat(bulks).to_csv(O/'R3_historical_bulk_all_models.tsv',sep='\t',index=False)
(O/'LINEAGE_PAIRS_COMPLETE.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'raw_scRNA_samples':len(reg),'new_state_models':len(models),'paired_patients':int(pa.patient.nunique()),'old_results_preserved':True,'spec_sha256':hashlib.sha256((R/'REASSESSMENT_SPEC.json').read_bytes()).hexdigest()},indent=2));print('LINEAGE_PAIRS_COMPLETE',flush=True)
