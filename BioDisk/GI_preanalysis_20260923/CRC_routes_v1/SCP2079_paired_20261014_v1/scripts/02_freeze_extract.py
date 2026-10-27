from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,hashlib,datetime
from scipy import sparse
R=Path(__file__).resolve().parent.parent;O=R/'results';Q=R.parent/'SCP2079_qualification_20261014_v1'
qualification=json.loads((R/'MATRIX_QUALIFICATION.json').read_text(encoding='utf-8'))
assert qualification['scale'] in ['verified_raw_counts','ln1p_CP10k_author_normalized']
plan=json.loads((Q/'CONDITIONAL_ANALYSIS_PLAN.json').read_text(encoding='utf-8'))
scale=qualification['scale'];raw=scale=='verified_raw_counts'
spec={**plan,'status':'frozen_before_target_extraction','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matrix_qualification_sha256':hashlib.sha256((R/'MATRIX_QUALIFICATION.json').read_bytes()).hexdigest(),'scale':scale,'measurement_primary':'ln1p(patient-time pooled raw counts / all-gene pool *10000)' if raw else 'patient-time arithmetic mean of author ln1p_CP10k cell expression; no invented raw counts','measurement_sensitivity':None if raw else 'ln1p(patient-time mean of expm1(author ln1p_CP10k)); equal cell weighting, not UMI pooling','depth':'log1p(mean author total_counts per cell); original metadata UMI field, not nCount_RNA','populations':{'primary':'fixed_qc_ge20_all_sites','sensitivity_1':'fixed_qc_ge50_all_sites','sensitivity_2':'author_retained_ge20_all_sites','liver_descriptive':'fixed_qc_ge20_liver'},'models':{'M1':['depth','epithelial_identity','mature_colonocyte'],'M0':['depth']},'inferential_min_n':10,'inferential_min_residual_df':5,'prediction_min_training_residual_df':5,'HC3':'two-sided Student-t with model residual df','paired_change_test':'exact conditional sign-flip enumeration of the two-sided signed-rank statistic, dropping zero changes, with average ranks for ties','paired_change_family':'three readouts within population/measurement','family_for_target_slope':'two outcomes within population/measurement/model; retain unestimable hypotheses as p=1 only within correction calculation, not as reported results','minimum_gene_coverage':'all fixed panel genes required in this round; any missing/duplicated fixed symbol blocks extraction','target_expression_or_association_results_seen':False}
sp=R/'SPECIFICATION.json'
assert not sp.exists(),'Do not overwrite the frozen specification'
sp.write_text(json.dumps(spec,indent=2),encoding='utf-8')
(R/'SPEC_FREEZE.json').write_text(json.dumps({'sha256':hashlib.sha256(sp.read_bytes()).hexdigest(),'frozen_utc':spec['frozen_utc'],'target_results_seen':False},indent=2),encoding='utf-8')
X=sparse.load_npz(O/'matrix_verified.npz').tocsc()
feat=pd.read_csv(O/'matrix_features.tsv',sep='\t',header=None,dtype=str)
bar=pd.read_csv(O/'matrix_barcodes.tsv',sep='\t',header=None,dtype=str)[0]
m=pd.read_csv(Q/'metadata_all_cells.tsv',sep='\t').set_index('NAME').loc[bar].reset_index()
fixed=(m.total_counts>=500)&(m.nFeature_RNA>=200)&(m.pct_counts_mt<=25)&(m.doublet_score<=.3)
programs={**spec['controls'],**spec['outcomes']}
genes=sorted(set(['CGNL1']+sum(programs.values(),[])))
counts=feat[1].value_counts();assert all(counts.get(g,0)==1 for g in genes)
indices=[int(feat.index[feat[1].eq(g)][0]) for g in genes]
G=X[indices,:].tocsr()
rows=[]
for definition,mask in [('fixed_qc',fixed.to_numpy()),('author_retained',np.ones(len(m),dtype=bool))]:
 for (donor,time),ids in m.loc[mask].groupby(['donor_id','treatment']).groups.items():
  ids=np.array(list(ids));n=len(ids);part=G[:,ids]
  if raw:
   denom=float(X[:,ids].sum());v=np.log1p(np.asarray(part.sum(axis=1)).ravel()/denom*10000)
   measures=[('pooled_raw_CP10k',v)]
  else:
   v=np.asarray(part.mean(axis=1)).ravel()
   linear=part.copy();linear.data=np.expm1(linear.data)
   alt=np.log1p(np.asarray(linear.mean(axis=1)).ravel())
   measures=[('mean_author_log',v),('log_mean_author_linear',alt)]
  for measurement,v in measures:
   row={'definition':definition,'measurement':measurement,'donor_id':donor,'treatment':time,'tissueSite':m.loc[ids,'tissueSite'].iloc[0],'cells':n,'depth':float(np.log1p(m.loc[ids,'total_counts'].mean())),'total_UMI_metadata':float(m.loc[ids,'total_counts'].sum()),'mean_mito_percent':float(m.loc[ids,'pct_counts_mt'].mean())}
   row.update(dict(zip(genes,v)))
   for key,gg in programs.items():row[key]=float(np.mean([row[g] for g in gg]))
   rows.append(row)
profiles=pd.DataFrame(rows)
profiles.to_csv(O/'patient_time_profiles.tsv',sep='\t',index=False,encoding='utf-8')
D=[];E=[]
for (definition,measurement),p in profiles.groupby(['definition','measurement']):
 for donor,g in p.groupby('donor_id'):
  pre=g[g.treatment.eq('Pre')];on=g[g.treatment.eq('On')]
  E.append({'definition':definition,'measurement':measurement,'donor_id':donor,'tissueSite':g.tissueSite.iloc[0],'Pre_cells':0 if pre.empty else int(pre.cells.iloc[0]),'On_cells':0 if on.empty else int(on.cells.iloc[0]),'paired':len(pre)==len(on)==1})
  if len(pre)!=1 or len(on)!=1:continue
  a=pre.iloc[0];b=on.iloc[0];assert a.tissueSite==b.tissueSite
  row={'definition':definition,'measurement':measurement,'donor_id':donor,'tissueSite':a.tissueSite,'Pre_cells':int(a.cells),'On_cells':int(b.cells),'min_cells':min(int(a.cells),int(b.cells))}
  for key in ['depth','CGNL1']+list(programs):row[key]=float(b[key]-a[key])
  D.append(row)
pd.DataFrame(D).to_csv(O/'patient_deltas.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(E).to_csv(O/'patient_eligibility.tsv',sep='\t',index=False,encoding='utf-8')
(R/'EXTRACTION_COMPLETE.json').write_text(json.dumps({'scale':scale,'profile_rows':len(profiles),'delta_rows':len(D),'genes':genes,'target_zero_cells_retained':True,'originals_not_modified':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2),encoding='utf-8')
print('EXTRACTION_COMPLETE',len(profiles),len(D),scale,flush=True)