from pathlib import Path
import os,json,time,hashlib,datetime
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ.update(NO_PROXY='*',no_proxy='*',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8')
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import anndata as ad
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006')
B=ROOT.parent;BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
plan=json.loads((ROOT/'plan_v3.json').read_text());panels=plan['panels']
def read10x(p):
 with h5py.File(p) as h:
  g=h['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr()
  names=[v.decode() for v in g['features/name'][:]];barcodes=[v.decode() for v in g['barcodes'][:]]
 if len(set(names))<len(names):
  uniq=sorted(set(names));idx={s:i for i,s in enumerate(uniq)}
  X=X@sparse.csr_matrix((np.ones(len(names)),(range(len(names)),[idx[s] for s in names])),shape=(len(names),len(uniq)));names=uniq
 return X.astype(np.float32),np.array(names),barcodes
def scores(X,names):
 umi=np.asarray(X.sum(axis=1)).ravel();norm=sparse.diags(1e4/np.maximum(umi,1))@X;norm=norm.tocsr();norm.data=np.log1p(norm.data)
 out={};cov={}
 for k,gs in panels.items():
  ix=np.flatnonzero(np.isin(names,gs));cov[k]=len(ix)
  out[k]=np.asarray(norm[:,ix].mean(axis=1)).ravel() if len(ix)>=2 else np.full(X.shape[0],np.nan)
 return pd.DataFrame(out),umi,norm,cov
def standard(x):
 x=np.asarray(x,dtype=float);s=x.std()
 return (x-x.mean())/s if s>0 else np.zeros(len(x))
# Reuse frozen accepted-EC count caches for new panel measurement only.
reg=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text())
lib=[]
for r in reg:
 X=sparse.load_npz(r['path']).tocsr();gs=pd.read_csv(r['genes'],sep='\t')
 # Adapt explicit gene table column only; never select by target outcome.
 names=gs.iloc[:,0].astype(str).to_numpy()
 if X.shape[1]!=len(names) and X.shape[0]==len(names):X=X.T.tocsr()
 if X.shape[1]!=len(names):raise ValueError(('gene dimensions',r['sample'],X.shape,gs.shape))
 if len(set(names))<len(names):
  uniq=sorted(set(names));idx={s:i for i,s in enumerate(uniq)}
  X=X@sparse.csr_matrix((np.ones(len(names)),(range(len(names)),[idx[s] for s in names])),shape=(len(names),len(uniq)))
  names=np.array(uniq)
 total=np.asarray(X.sum(axis=0)).ravel();n=X.shape[0]
 profile=np.log1p(total/max(total.sum(),1)*1e4)
 z={'dataset':r['dataset'],'sample':r['sample'],'patient':r['donor'],'n_accepted_EC':n,'total_UMI':float(total.sum()),'eligible_v3_n20':bool(n>=20)}
 for k,panel in {'CGNL1':['CGNL1'],**panels}.items():
  ix=np.flatnonzero(np.isin(names,panel));z[k]=float(profile[ix].mean()) if len(ix)>=(1 if k=='CGNL1' else 2) else None;z[k+'_covered']=len(ix)
 lib.append(z)
libs=pd.DataFrame(lib);libs.to_csv(ROOT/'results/accepted_EC_library_measurements.tsv',sep='\t',index=False)
elig=libs[libs.eligible_v3_n20].copy()
if len(elig):
 elig.groupby(['dataset','patient'],as_index=False)[['CGNL1','notch_response','adhesion','hypoxia','historical_junction']].mean().to_csv(ROOT/'results/accepted_EC_patient_measurements.tsv',sep='\t',index=False)

print('GENE_COLLAPSE_REPAIRED')
