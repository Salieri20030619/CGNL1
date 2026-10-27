from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
import json,gzip,io,itertools,datetime,traceback
import pandas as pd,numpy as np,h5py
from scipy import sparse,io as sio,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import anndata as ad
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'feasibility_audit_20261008_v1';O=R/'results';BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
checks=[]
def check(name,ok,detail):
 checks.append(dict(check=name,passed=bool(ok),detail=detail));print(name,bool(ok),flush=True)
def table(name,rows):
 t=pd.DataFrame(rows);t.to_csv(O/(name+'.tsv'),sep='\t',index=False);return t
rows=[]
for p in sorted((B/'spatial_v3/raw').glob('GSM981545*_cell_feature_matrix.h5')):
 with h5py.File(p) as h:
  g=h['matrix'];types=np.array([v.decode() for v in g['features/feature_type'][:]]);names=np.array([v.decode() for v in g['features/name'][:]]);ids=np.array([v.decode() for v in g['features/id'][:]])
  target=np.isin(names,['CGNL1','JACOP','PCING','ENSG00000128849','84952'])|np.isin(ids,['CGNL1','JACOP','PCING','ENSG00000128849','84952'])
  rr=dict(file=p.name,features=len(names),gene_expression_features=int((types=='Gene Expression').sum()),CGNL1_measured=bool(target.any()),raw_cells=int(g['shape'][1]))
  rows.append(rr)
  if p.name.startswith('GSM9815454'):
   X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr();bc=[x.decode() for x in g['barcodes'][:]]
   with gzip.open(BIG/'downloads/GSM9815454_cells.parquet.gz','rb') as f:obs=pd.read_parquet(io.BytesIO(f.read())).set_index('cell_id').loc[bc]
   gene=np.asarray(X[:,types=='Gene Expression'].sum(1)).ravel();allc=np.asarray(X.sum(1)).ravel()
   check('M242T_raw_counts_metadata_alignment_corrected_validator',np.array_equal(gene,obs.transcript_counts),{'gene_count_max_diff':float(np.abs(gene-obs.transcript_counts).max()),'all_feature_vs_transcript_mismatched_cells':int((allc!=obs.transcript_counts).sum()),'reason_first_audit_false':'The first audit compared all features including negative controls with gene-only transcript_counts; original analysis correctly selects Gene Expression. This is an audit-validator error, not an original data error.'})
table('all_four_Xenium_measurement_audit',rows)

# Source-to-cache check on first nonempty sample by accession in each study, fixed without target selection.
reg=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text());rawchecks=[]
for ds in sorted({x['dataset'] for x in reg}):
 rec=sorted([x for x in reg if x['dataset']==ds and x['selected_cells']>0],key=lambda r:r['sample'])[0]
 p=Path(rec['source_matrix']);stem=str(p).replace('matrix.mtx.gz','')
 feats=pd.read_csv(stem+'features.tsv.gz',sep='\t',header=None);bc=pd.read_csv(stem+'barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
 meta=pd.read_csv(rec['cells'],sep='\t');indices=pd.Index(bc).get_indexer(meta.barcode)
 print('SOURCE_READ',rec['sample'],flush=True)
 with gzip.open(p,'rb') as f,threadpool_limits(limits=4):X=sio.mmread(f).tocsr()
 Y=sparse.load_npz(rec['path']).tocsr();diff=X[:,indices]-Y
 ok=(indices>=0).all() and diff.nnz==0 and np.array_equal(feats[1].astype(str),pd.read_csv(rec['genes'],sep='\t').iloc[:,0].astype(str))
 rawchecks.append(dict(dataset=ds,sample=rec['sample'],selected_cells=len(indices),raw_shape=str(X.shape),different_entries=diff.nnz,passed=bool(ok)))
 check(rec['sample']+'_raw_MTX_to_cache',ok,rawchecks[-1]);del X,Y
table('raw_MTX_cache_spotchecks',rawchecks)

ar=B/'A04_subtype_check_20261002_v1/results'
pt=pd.read_csv(ar/'patient_profiles_composition.tsv',sep='\t');old=pd.read_csv(ar/'composition_adjusted_models.tsv',sep='\t');reps=[]
for _,r in old.iterrows():
 if r.status!='estimated':continue
 v=pt[pt.branch==r.branch].copy();F=pd.DataFrame(index=v.index)
 for col in r.covariates.split('|'):
  F[col]=1. if col=='const' else (v.dataset==col).astype(float) if col.startswith('GSE') else v[col]
 m=sm.OLS(v.junction,F).fit(use_t=True);beta=m.params.CGNL1*v.CGNL1.std()/v.junction.std()
 reps.append(dict(branch=r.branch,model=r.model,family=r.family,n=len(v),beta=beta,P=m.pvalues.CGNL1,old_FDR=r.FDR,pass_beta=abs(beta-r.beta)<1e-8,pass_P=abs(m.pvalues.CGNL1-r.P)<1e-8))
t=pd.DataFrame(reps);t['FDR']=np.nan
for fam,ix in t.groupby('family').groups.items():t.loc[ix,'FDR']=multipletests(t.loc[ix,'P'],method='fdr_bh')[1]
t['pass_FDR']=abs(t.FDR-t.old_FDR)<1e-8;table('A04_composition_reproduction',t);check('A04_all_composition_models_reproduced',t[['pass_beta','pass_P','pass_FDR']].to_numpy().all(),{'models':len(t)})

for prefix,path,modelfile,patientfile,keys in [
 ('A06',B/'review_consensus_retest_20261002_v1/results','matched_models.tsv','matched_patient_effects.tsv',['stage']),
 ('A04_same_subtype',ar,'same_subtype_matched_models.tsv','matched_patient_effects.tsv',['branch','constraint','minpairs'])]:
 old=pd.read_csv(path/modelfile,sep='\t');pt=pd.read_csv(path/patientfile,sep='\t');rows=[]
 for _,r in old.iterrows():
  mask=np.ones(len(pt),bool)
  for k in keys:mask &= pt[k].eq(r[k]).to_numpy()
  if r.scope!='pooled':mask &= pt.dataset.eq(r.scope).to_numpy()
  x=pt.loc[mask,'effect'].to_numpy();pv=float(stats.ttest_1samp(x,0).pvalue) if len(x)>=4 and np.std(x)>0 else np.nan
  eff=float(x.mean()) if len(x) else np.nan
  rows.append({**{k:r[k] for k in keys},'scope':r.scope,'family':r.family,'n':len(x),'effect':eff,'P':pv,'old_FDR':r.FDR,'pass_effect':bool(np.isclose(eff,r.effect,equal_nan=True)),'pass_P':bool(np.isclose(pv,r.P,equal_nan=True))})
 t=pd.DataFrame(rows);t['FDR']=np.nan
 for fam,ix in t.groupby('family').groups.items():
  valid=t.loc[ix].index[t.loc[ix,'P'].notna()]
  if len(valid):t.loc[valid,'FDR']=multipletests(t.loc[valid,'P'],method='fdr_bh')[1]
 t['pass_FDR']=np.isclose(t.FDR,t.old_FDR,equal_nan=True);table(prefix+'_matched_model_reproduction',t);check(prefix+'_matched_models_reproduced',t[['pass_effect','pass_P','pass_FDR']].to_numpy().all(),{'models':len(t)})

# Re-evaluate native LIANA at the corrected physical scale, retaining all original pairs and expression thresholds.
import scanpy as sc,liana as li
cfg=json.loads((B/'v3_pilot_20261006/results/LIANA_PILOT_PARAMETERS.json').read_text());resource=pd.DataFrame(cfg['pairs']);outrows=[]
for rec in json.loads((B/'v3_pilot_20261006/results/TECHNICAL_INPUTS.json').read_text()):
 a=ad.read_h5ad(rec['path']);a.obsm['spatial']=a.obsm['spatial']*65/55;sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 li.utils.spatial_neighbors(a,bandwidth=200,kernel='linear',cutoff=0,max_neighbours=100,standardize=True,key_added='audit_liana')
 key=[k for k in a.obsp if k.startswith('audit_liana')][0]
 result=li.mt.bivariate(a,resource=resource,local_name='cosine',global_name=None,connectivity_key=key,n_perms=None,nz_prop=.01,use_raw=False,verbose=False)
 dense=result.X.toarray() if sparse.issparse(result.X) else np.asarray(result.X)
 old=pd.read_csv(B/'v3_pilot_20261006/results'/f'{rec["sample"]}_LIANA_pair_scores.tsv',sep='\t').set_index('pair')
 original=ad.read_h5ad(BIG/'inputs'/f'{rec["sample"]}_LIANA_scores.h5ad');original=original[:,result.var_names];before=original.X.toarray() if sparse.issparse(original.X) else np.asarray(original.X)
 for i,pair in enumerate(result.var_names):outrows.append(dict(sample=rec['sample'],pair=pair,old_mean=float(old.loc[pair,'mean_score']),corrected_mean=float(dense[:,i].mean()),location_score_correlation=float(np.corrcoef(before[:,i],dense[:,i])[0,1]),mean_absolute_change=float(np.abs(before[:,i]-dense[:,i]).mean()),P=np.nan,FDR=np.nan,scope='technical_scores_only_no_cell_specific_causality'))
 check(rec['sample']+'_corrected_LIANA_native_run',np.isfinite(dense).all(),{'pairs':result.n_vars,'coordinate_multiplier':65/55})
table('LIANA_corrected_scale_retest',outrows)
(O/'FOLLOWUP_CHECKS.json').write_text(json.dumps(checks,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
(O/'FOLLOWUP_COMPLETE.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks),'superseded_check':'M242T_raw_counts_metadata_alignment: validator must use Gene Expression counts'},indent=2))
print('FOLLOWUP_COMPLETE',flush=True)
