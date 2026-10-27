from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,gzip,hashlib,datetime,sys
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy import sparse
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent.parent
Q=R.parent/'SCP2079_qualification_20261014_v1'
O=R/'results';O.mkdir(exist_ok=True)
def save(name,obj): (R/name).write_text(json.dumps(obj,indent=2),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
P=R/'raw/expression_epi_counts'
meta=pd.read_csv(Q/'metadata_all_cells.tsv',sep='\t',dtype={'NAME':str,'donor_id':str,'sampleID':str,'treatment':str})
bar=pd.read_csv(P/'barcodes.tsv.gz',sep='\t',header=None,dtype=str)[0]
feat=pd.read_csv(P/'features.tsv.gz',sep='\t',header=None,dtype=str)
assert meta.NAME.is_unique and bar.is_unique
assert set(bar)==set(meta.NAME)
meta=meta.set_index('NAME').loc[bar].reset_index()
assert feat[0].is_unique
with gzip.open(P/'matrix.mtx.gz','rt',encoding='utf-8') as f:header=[next(f).strip() for _ in range(5)]
print('HEADER',header,flush=True)
print('READING_MATRIX',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
with threadpool_limits(limits=4):
 with gzip.open(P/'matrix.mtx.gz','rb') as f: coo=mmread(f)
print('READ_COMPLETE',coo.shape,coo.nnz,flush=True)
assert coo.shape==(len(feat),len(bar))
nnz_input=coo.nnz
X=coo.tocsc();del coo
assert X.nnz==nnz_input,'Duplicated coordinates'
negative=nonfinite=noninteger=zeros=0
for start in range(0,len(X.data),5000000):
 a=X.data[start:start+5000000]
 negative+=int((a<0).sum());nonfinite+=int((~np.isfinite(a)).sum());zeros+=int((a==0).sum());noninteger+=int((np.abs(a-np.rint(a))>1e-7).sum())
assert negative==nonfinite==0
sums=np.asarray(X.sum(axis=0)).ravel()
def match(a,b):
 d=np.abs(a-b)
 return {'max_abs':float(d.max()),'max_relative':float((d/np.maximum(1,np.abs(b))).max()),'median_relative':float(np.median(d/np.maximum(1,np.abs(b)))),'fraction_close':float(np.isclose(a,b,rtol=1e-5,atol=1e-5).mean())}
vs_total=match(sums,meta.total_counts.to_numpy())
vs_ncount=match(sums,meta.nCount_RNA.to_numpy())
scale='verified_raw_counts' if noninteger==0 and vs_total['fraction_close']>.99999 else 'unresolved'
exp_sum=None
if noninteger:
 # Test numerical units without reconstructing or rounding raw counts.
 exp_sum=np.zeros(X.shape[1])
 for j in range(X.shape[1]):
  a=X.data[X.indptr[j]:X.indptr[j+1]]
  exp_sum[j]=np.expm1(a).sum()
 if vs_ncount['fraction_close']>.999 and np.max(np.abs(exp_sum-10000))<2:scale='ln1p_CP10k_author_normalized'
actual_features=np.diff(X.indptr)
feature_match=match(actual_features,meta.nFeature_RNA.to_numpy())
qc=(meta.total_counts>=500)&(meta.nFeature_RNA>=200)&(meta.pct_counts_mt<=25)&(meta.doublet_score<=.3)
profiles=meta[['NAME','donor_id','treatment','tissueSite','sampleID','total_counts','nCount_RNA','nFeature_RNA','pct_counts_mt','doublet_score']].copy()
profiles['matrix_sum']=sums;profiles['matrix_detected_genes']=actual_features;profiles['fixed_qc']=qc
if exp_sum is not None:profiles['expm1_matrix_sum']=exp_sum
profiles.to_csv(O/'matrix_cell_qualification.tsv',sep='\t',index=False,encoding='utf-8')
summary={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matrix_path':str(P/'matrix.mtx.gz'),'shape':list(X.shape),'nnz':int(X.nnz),'noninteger_nonzeros':noninteger,'negative':negative,'nonfinite':nonfinite,'explicit_zeros':zeros,'cell_barcodes_match_full_metadata':True,'features_unique':True,'symbol_duplicates':int(feat[1].duplicated().sum()),'matrix_sum_vs_total_counts':vs_total,'matrix_sum_vs_nCount_RNA':vs_ncount,'detected_genes_vs_nFeature_RNA':feature_match,'scale':scale,'expm1_column_sum_quantiles':None if exp_sum is None else np.quantile(exp_sum,[0,.25,.5,.75,1]).tolist(),'header':header,'files':[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(P.iterdir())]}
save('MATRIX_QUALIFICATION.json',summary)
# Save a server-only sparse cache once, never on the local computer.
sparse.save_npz(O/'matrix_verified.npz',X,compressed=False)
feat.to_csv(O/'matrix_features.tsv',sep='\t',header=False,index=False,encoding='utf-8')
bar.to_csv(O/'matrix_barcodes.tsv',sep='\t',header=False,index=False,encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
assert scale!='unresolved','Scale unresolved: stop before target extraction or inference'