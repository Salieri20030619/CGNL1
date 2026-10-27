from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,hashlib,gzip,tarfile,io
import numpy as np,pandas as pd
from scipy import io as sio,sparse
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from threadpoolctl import threadpool_limits
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');O=Path.cwd()
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
panels=json.loads((R/'epithelial_independence_20261012_v1/MEASUREMENT_AMENDMENT01.json').read_text())['common_measured_panels'];panels['mural']=['RGS5','CSPG4','MCAM','PDGFRB','NOTCH3']
comparators=['CDH1','KRT20','CLDN3','CLDN4','VIL1','CDX2'];targets=['CGNL1']+comparators
spec={'frozen_utc':now(),'question':'Does CGNL1 expression in pathology-defined macrometastatic tumor regions associate with proximity to hepatic tissue and with interface-specific mural neighborhood organization beyond fixed local-state and lineage proxies?',
 'estimand_limit':'regional mixed-spot marker, not cell-intrinsic expression, true physical contact, HGP label or causal mechanism',
 'prior_knowledge':'old unadjusted distance, mural-neighborhood and micro/macro results known; new restricted interface estimands not yet calculated; exploratory reuse of cohort',
 'population':'GSE294385 non-recurrence liver sections; per-section author Layer2/Layer3-consistent macrometastatic tumor, UMI>=500 and >=200 genes',
 'distance':'nearest author Hepatic lobule spot using coordinates on the measured tissue grid; target-excluded labels; proximity proxy, not exact histological interface',
 'geometry':'macro components at adjacent 110um; each retained component must contain >=5 near and >=5 deep macro spots; every eligible section >=20 near and >=20 deep; compare within section, aggregate equal section weights within patient',
 'primary_bands':{'near':[50,200],'deep':[500,1500]},'sensitivity_bands':{'near':[100,300],'deep':[600,1500]},
 'panels':panels,'targets':targets,'main_neighbor_radius_um':200,'radius_sensitivity_um':[100,300],
 'main_covariates':['logUMI','epithelial_identity','mature_colonocyte','fibroblast','endothelial','hepatocyte','mural','hypoxia'],
 'E1':'CGNL1 ln1pCP10k near-minus-deep coefficient conditioned on fixed covariates and macro-component intercepts',
 'E2':'CGNL1 association with 200um neighbor mural score differs between near and deep; coefficient of near x within-section standardized neighbor_mural, conditional on same covariates, near and neighbor_mural main effects and component intercepts',
 'secondary_model':'E1 adjusted for depth only and component; descriptive decomposition, not replacement if main fails',
 'specificity':'repeat both endpoints for six fixed ordinary epithelial genes; compare CGNL1 within-section-SD effect to median comparator effect at patient level',
 'composition_matching':'sensitivity match near/deep without target values on main covariates within section/component with standardized RMS distance <=0.75, without replacement; report retained patients and support; do not force estimate outside overlap',
 'inference':'patient equal-weight means; two-sided exact sign-flip test on patient effects; n<10 exploratory with valid independent patient repetition; n<2 no P; no spot P; t CI descriptive with small-n caveat; leave-one-patient effect range',
 'FDR':'BH on two CGNL1 main endpoint hypotheses; failed hypothesis retained as P=1 for family size only, actual P missing; sensitivity branch families separate and labeled; comparator family separate; no selection by FDR',
 'decision':'nominal raw P<0.05 with effect/concordance and robustness; novelty and independent source qualification also required for advancing a definitive mainline; negative permitted',
 'spatial_validation':'patient aggregation avoids treating spots as independent; complementary spatial-block aggregation and local-neighborhood conditional sensitivity; no claim that conventional random spot permutation is a valid null',
 'no_new_cohort_identity_claim':'GSE217414 lacks same pathology masks, so existing external generic mural association cannot validate E1/E2',
 'primary_effects_seen_before_freeze':False}
f=O/'ANALYSIS_SPECIFICATION.json';assert not f.exists();b=(json.dumps(spec,ensure_ascii=False,indent=2)+'\n').encode('utf-8');f.write_bytes(b);(O/'SPEC_FREEZE.json').write_text(json.dumps({'utc':now(),'sha256':hashlib.sha256(b).hexdigest()},indent=2),encoding='utf-8')
ann=pd.read_csv(R/'results/spatial_author_metadata.tsv.gz',sep='\t');ann['sample']=ann.Barcode.str.split('~').str[0];ann['barcode']=ann.Barcode.str.split('~').str[1]
legacy=pd.read_csv(R/'spatial_v3/results/screen/spot_target_coordinates.tsv.gz',sep='\t');jobs=[j for j in json.loads((R/'spatial_v3/metadata/LIVER_DOWNLOAD_PLAN.json').read_text()) if 'Recurrence' not in j['sample_type']]
wanted=sorted(set(targets+sum(panels.values(),[])));tables=[];coverage=[];manifest=[]
for j in jobs:
 name=j['sample'];path=R/'spatial_v3/raw/liver'/f'{name}.tar.gz';print('EXTRACT',name,flush=True)
 with tarfile.open(path,'r|gz') as tar:
  for m in tar:
   if not m.isfile():continue
   if m.name.endswith('filtered_feature_bc_matrix/features.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as ff:genes=pd.read_csv(ff,sep='\t',header=None)[1].astype(str).to_numpy()
   elif m.name.endswith('filtered_feature_bc_matrix/barcodes.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as ff:bc=pd.read_csv(ff,sep='\t',header=None)[0].astype(str).to_numpy()
   elif m.name.endswith('filtered_feature_bc_matrix/matrix.mtx.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as ff,threadpool_limits(limits=4):X=sio.mmread(ff).tocsr()
 a=ann[ann['sample'].eq(name)].copy().reset_index(drop=True);ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all();X=X[:,ix];X.eliminate_zeros();umi=np.asarray(X.sum(0)).ravel();ng=X.getnnz(0)
 old=legacy[legacy['sample'].eq(name)].set_index('barcode').loc[a.barcode];xy=old[['x_um','y_um']].to_numpy();assert np.array_equal(umi,old.UMI)
 vals={};counts={}
 for gene in wanted:
  indices=np.flatnonzero(genes==gene)
  if len(indices):cnt=np.asarray(X[indices].sum(0)).ravel();counts[gene]=cnt;vals[gene]=np.log1p(cnt/np.maximum(umi,1)*1e4)
 assert np.array_equal(counts['CGNL1'],old.CGNL1_counts)
 good=(umi>=500)&(ng>=200);idx=np.flatnonzero(good);v=pd.DataFrame({'sample':name,'patient':j['patient'],'barcode':a.barcode.to_numpy()[idx],'Layer2':a.Layer2.to_numpy()[idx],'region':a.Layer3.to_numpy()[idx],'x_um':xy[idx,0],'y_um':xy[idx,1],'UMI':umi[idx],'logUMI':np.log1p(umi[idx]),'n_genes':ng[idx]})
 hep=a.Layer3.eq('Hepatic lobule').to_numpy();dist=cKDTree(xy[hep]).query(xy[idx])[0] if hep.any() else np.full(len(idx),np.nan);v['distance']=dist
 for gene in targets:v[gene]=vals[gene][idx] if gene in vals else np.nan
 v['CGNL1_count']=counts['CGNL1'][idx]
 for panel,gs in panels.items():
  available=[g for g in gs if g in vals];coverage.append({'sample':name,'panel':panel,'measured':len(available),'expected':len(gs),'missing':'|'.join(g for g in gs if g not in vals)})
  v[panel]=np.mean([vals[g][idx] for g in available],axis=0) if len(available)>=.8*len(gs) else np.nan
 tree=cKDTree(xy[idx]); pair=tree.query_pairs(110+1e-6,output_type='ndarray')
 macro=v.region.eq('Liver macrometastasis tumor').to_numpy();pairs=pair[macro[pair[:,0]]&macro[pair[:,1]]];G=sparse.csr_matrix((np.ones(len(pairs)*2),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(v),len(v)));_,labels=connected_components(G,directed=False);v['macro_component']=np.where(macro,labels,-1)
 for radius in [100,200,300]:
  pairs=tree.query_pairs(radius+1e-6,output_type='ndarray');rr=np.r_[pairs[:,0],pairs[:,1]];cc=np.r_[pairs[:,1],pairs[:,0]];G=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(v),len(v)));degree=np.asarray(G.sum(1)).ravel();W=sparse.diags(1/np.maximum(degree,1))@G
  v['degree_'+str(radius)]=degree
  v['neighbor_mural_'+str(radius)]=np.where(degree>0,W@v.mural.to_numpy(),np.nan)
  for panel in ['fibroblast','endothelial','hepatocyte','epithelial_identity']:
   if radius==200:v['neighbor_'+panel]=np.where(degree>0,W@v[panel].to_numpy(),np.nan)
 tables.append(v)
 manifest.append({'sample':name,'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'QC_spots':len(v),'matrix_raw_UMI_matches_previous':True,'CGNL1_raw_counts_match_previous':True})
 del X
d=pd.concat(tables,ignore_index=True);d.to_csv(O/'results/spatial_spot_features.tsv.gz',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'metadata/panel_coverage.tsv',sep='\t',index=False)
(O/'metadata/RAW_SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
summary={'utc':now(),'sections':len(tables),'patients':int(d.patient.nunique()),'QC_spots':len(d),'model_effects_not_yet_calculated':True,'layers':d[['Layer2','region']].drop_duplicates().to_dict('records')};(O/'EXTRACTION_COMPLETE.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary),flush=True)
