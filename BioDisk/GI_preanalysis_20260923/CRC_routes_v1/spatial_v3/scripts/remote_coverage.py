from pathlib import Path
import sys,json,h5py,numpy as np,pandas as pd
from scipy.sparse import csc_matrix
sys.stdout.reconfigure(encoding='utf-8')
O=Path(sys.argv[1]);rows=[]
for p in sorted((O/'raw').glob('*.h5')):
 with h5py.File(p) as f:
  m=f['matrix'];n=[x.decode() if isinstance(x,bytes) else str(x) for x in m['features/name'][:]]
  ids=[x.decode() if isinstance(x,bytes) else str(x) for x in m['features/id'][:]]
  shape=m['shape'][:];ix=[i for i,x in enumerate(n) if x=='CGNL1' or ids[i].split('.')[0]=='ENSG00000128849'];d={'file':p.name,'features':len(n),'spots_or_cells':int(shape[1]),'CGNL1_measured':bool(ix)}
  if ix:
   X=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=shape);v=np.asarray(X[ix,:].sum(axis=0)).ravel();d.update(CGNL1_total_counts=int(v.sum()),CGNL1_positive_spots=int((v>0).sum()),CGNL1_detected_pct=float((v>0).mean()*100))
  rows.append(d)
pd.DataFrame(rows).to_csv(O/'results/CGNL1_platform_coverage.tsv',sep='\t',index=False);print(pd.DataFrame(rows).to_string(index=False))