from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import json,datetime
import numpy as np,pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');O=Path.cwd();d=pd.read_csv(R/'spatial_v3/results/screen/spot_target_coordinates.tsv.gz',sep='\t');e=pd.read_csv(O/'results/model_eligibility.tsv',sep='\t');samples=e[(e.variant=='primary')&(e.endpoint=='E1')&(e.status=='estimated')]['sample'].unique();out=[]
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'qualify nearest hepatic spot as interface proxy, after target results; do not change primary definition','checks':['Layer2/Layer3 consistent hepatic anchors','geodesic distance on author tissue grid using adjacent110um edges, no crossing unmapped holes'],'effects':'no new expression effect examined; check band assignment changes only'};(O/'GEOMETRY_DIAGNOSTIC_SPEC.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
for sample in samples:
 v=d[d['sample']==sample].reset_index(drop=True);xy=v[['x_um','y_um']].to_numpy();macro=v.Layer3.eq('Liver macrometastasis tumor')&v.Layer2.eq('Liver_metastatic carcinoma')&v.QC_pass.astype(str).str.lower().isin(['true','1']);hep=v.Layer3.eq('Hepatic lobule')&v.Layer2.eq('Non-neoplastic liver');old=v.distance_to_hepatic_um.to_numpy();eu=cKDTree(xy[hep]).query(xy)[0];pairs=cKDTree(xy).query_pairs(110+1e-6,output_type='ndarray');dist=np.sqrt(((xy[pairs[:,0]]-xy[pairs[:,1]])**2).sum(1));G=sparse.csr_matrix((np.r_[dist,dist],(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(v),len(v)));geo=dijkstra(G,directed=False,indices=np.flatnonzero(hep),min_only=True)
 def band(x):return np.where((x>=50)&(x<=200),1,np.where((x>=500)&(x<=1500),2,0))
 a=band(old);b=band(eu);c=band(geo);sel=macro.to_numpy();ret=sel&(a>0)
 out.append(dict(sample=sample,patient=int(v.patient.iloc[0]),macro_QC=int(sel.sum()),primary_band_spots=int(ret.sum()),anchor_consistency_changed=int(((a!=b)&ret).sum()),anchor_consistency_changed_fraction=float(((a!=b)&ret).sum()/ret.sum()),geodesic_changed=int(((a!=c)&ret).sum()),geodesic_changed_fraction=float(((a!=c)&ret).sum()/ret.sum()),disconnected_macro=int((sel&~np.isfinite(geo)).sum()),geo_near=int((sel&(c==1)).sum()),geo_deep=int((sel&(c==2)).sum())))
q=pd.DataFrame(out);q.to_csv(O/'results/geometry_diagnostic.tsv',sep='\t',index=False);print(q.to_string(index=False));(O/'GEOMETRY_DIAGNOSTIC_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'max_anchor_change':float(q.anchor_consistency_changed_fraction.max()),'max_geodesic_change':float(q.geodesic_changed_fraction.max())},indent=2),encoding='utf-8')
