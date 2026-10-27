from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import ast,json,hashlib,datetime,itertools
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');O=Path.cwd();S=json.loads((O/'ANALYSIS_SPECIFICATION.json').read_text());targets=S['targets'];covs=S['main_covariates']
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'post-primary geometry sensitivity','reason':'geodesic diagnostic changed band labels in up to31.5% of primary-band macro spots; one section has disconnected tissue','change':'replace Euclidean nearest-hepatic distance by shortest path on 110um-adjacent author-tissue grid; anchor Hepatic lobule and Non-neoplastic liver consistent labels; unreachable spots remain ineligible','unchanged':'target, comparator panel, covariates, near50-200/deep500-1500 bands, component rules, two endpoints and patient-unit tests','not_new_independent_validation':True}
(O/'GEODESIC_VALIDATION_SPEC.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
source=(O/'scripts/cgnl_05_spatial_models.py').read_text(encoding='utf-8');tree=ast.parse(source);functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['select','design','infer']];exec(compile(ast.Module(body=functions,type_ignores=[]),'<frozen model functions>','exec'))
d=pd.read_csv(O/'results/spatial_spot_features.tsv.gz',sep='\t');old=pd.read_csv(R/'spatial_v3/results/screen/spot_target_coordinates.tsv.gz',sep='\t');rows=[];eligible=[];qa=[]
for sample,v in d.groupby('sample'):
 a=old[old['sample']==sample].reset_index(drop=True);xy=a[['x_um','y_um']].to_numpy();hep=a.Layer3.eq('Hepatic lobule')&a.Layer2.eq('Non-neoplastic liver');pairs=cKDTree(xy).query_pairs(110+1e-6,output_type='ndarray');dd=np.sqrt(((xy[pairs[:,0]]-xy[pairs[:,1]])**2).sum(1));G=sparse.csr_matrix((np.r_[dd,dd],(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(a),len(a)));dist=dijkstra(G,directed=False,indices=np.flatnonzero(hep),min_only=True);mapping=pd.Series(dist,index=a.barcode);v=v.copy();v['distance']=mapping.loc[v.barcode].to_numpy();v=select(v,50,200,500,1500).dropna(subset=covs+targets+['neighbor_mural_200']);patient=int(d[d['sample']==sample].patient.iloc[0]);info=dict(sample=sample,patient=patient,near=int(v.near.sum()),deep=int(len(v)-v.near.sum()))
 if min(info['near'],info['deep'])<20:eligible.append(dict(info,status='insufficient_band_support'));continue
 for endpoint in ['E1','E2']:
  X,index,names=design(v,covs,200,endpoint)
  if X is None or np.linalg.matrix_rank(X)<X.shape[1] or len(X)-X.shape[1]<10:eligible.append(dict(info,endpoint=endpoint,status='not_estimable'));continue
  Y=v[targets].to_numpy();b=np.linalg.lstsq(X,Y,rcond=None)[0][index];N=np.delete(X,index,axis=1);rx=X[:,index]-N@np.linalg.lstsq(N,X[:,index],rcond=None)[0];ry=Y-N@np.linalg.lstsq(N,Y,rcond=None)[0];qa.append(float(np.max(abs(b-rx@ry/(rx@rx)))));eligible.append(dict(info,endpoint=endpoint,status='estimated'))
  for g,x,ys in zip(targets,b,Y.std(0)):rows.append(dict(sample=sample,patient=patient,endpoint=endpoint,gene=g,effect=float(x),effect_SD=float(x/ys)))
q=pd.DataFrame(rows);q.to_csv(O/'results/geodesic_section_effects.tsv',sep='\t',index=False);pd.DataFrame(eligible).to_csv(O/'results/geodesic_eligibility.tsv',sep='\t',index=False);p=q.groupby(['endpoint','gene','patient'],as_index=False)[['effect','effect_SD']].mean();p.to_csv(O/'results/geodesic_patient_effects.tsv',sep='\t',index=False);summ=[]
for (endpoint,gene),g in p.groupby(['endpoint','gene']):summ.append(dict(endpoint=endpoint,gene=gene,**infer(g.effect),mean_SD=float(g.effect_SD.mean())))
t=pd.DataFrame(summ);t['FDR']=np.nan
for role,idx in t.assign(role=np.where(t.gene=='CGNL1','target','comparators')).groupby('role').groups.items():t.loc[idx,'FDR']=multipletests(t.loc[idx,'P'].fillna(1),method='fdr_bh')[1]
t.to_csv(O/'results/geodesic_summary.tsv',sep='\t',index=False);specific=[]
for endpoint,g in p.groupby('endpoint'):
 mat=g.pivot(index='patient',columns='gene',values='effect_SD').reindex(columns=targets).dropna();specific.append(dict(endpoint=endpoint,**infer(mat.CGNL1-mat[targets[1:]].median(axis=1))))
sp=pd.DataFrame(specific);sp['FDR']=multipletests(sp.P.fillna(1),method='fdr_bh')[1];sp.to_csv(O/'results/geodesic_specificity.tsv',sep='\t',index=False);assert max(qa)<1e-8
(O/'GEODESIC_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sections':q['sample'].nunique(),'patients':p.patient.nunique(),'FWL_max_error':max(qa)},indent=2),encoding='utf-8');print('TARGET',t[t.gene=='CGNL1'].to_string(index=False));print('SPECIFICITY',sp.to_string(index=False))
