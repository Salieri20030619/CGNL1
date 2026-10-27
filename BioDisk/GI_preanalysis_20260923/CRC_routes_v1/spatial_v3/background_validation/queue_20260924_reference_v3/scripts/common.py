import pathlib,json,sys,gzip,hashlib
import numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
Q=pathlib.Path(sys.argv[1]);C=json.loads((Q/'CONFIG.json').read_text());I=pathlib.Path(C['input_run']);S=pathlib.Path(C['spatial_root']);DEEP=S.parent/'deep_v2'
def dump(path,obj):
 path=pathlib.Path(path);path.parent.mkdir(exist_ok=True,parents=True);t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()),encoding='utf-8');t.replace(path)
def load_sample(sample):
 get=lambda suffix:I/next(p for p in sample['files'] if p.endswith(suffix));ft=pd.read_csv(get('_features.tsv.gz'),sep='\t',header=None);bc=pd.read_csv(get('_barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str).to_numpy()
 with gzip.open(get('_matrix.mtx.gz'),'rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsr()
 X.eliminate_zeros();assert X.shape==(len(ft),len(bc));assert pd.Index(bc).is_unique;lib=np.asarray(X.sum(axis=0)).ravel();ng=np.asarray(X.getnnz(axis=0)).ravel();sy=ft[1].astype(str).to_numpy();mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;good=(lib>=500)&(ng>=200)&(mt<=25);return X[:,good],sy,bc[good],lib[good]
def role(s):
 if s['dataset']=='GSE285991':return 'liver_tumor'
 if s['title'].startswith('L'):return 'liver_normal' if s['title'].endswith('N') else 'liver_tumor'
 return 'colon_normal' if s['title'].endswith('N') else 'colon_tumor'
def bh(p,n=None):
 p=np.asarray(p,float);orig=len(p);n=n or orig;p=np.r_[np.nan_to_num(p,nan=1),np.ones(n-orig)];o=np.argsort(p);q=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[o]=np.minimum(q,1);return out[:orig]
