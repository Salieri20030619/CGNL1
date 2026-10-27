from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4')
import json, hashlib, datetime, traceback, gzip, io
import numpy as np
import pandas as pd
import h5py, anndata as ad
from scipy import sparse, stats
from scipy.spatial import cKDTree, distance
from sklearn.metrics import adjusted_rand_score
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'feasibility_audit_20261008_v1'
P=B/'v3_pilot_20261006'
BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
R.mkdir(exist_ok=True)
OUT=R/'results'; OUT.mkdir(exist_ok=True)
plan=json.loads((P/'plan_v3.json').read_text()); panels=plan['panels']
checks=[]; fingerprints={}; modules={}
def sha(p):
 p=Path(p); h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''): h.update(chunk)
 fingerprints[str(p)]={'sha256':h.hexdigest(),'bytes':p.stat().st_size}
 return h.hexdigest()
def chk(name,passed,detail):
 checks.append(dict(check=name,passed=bool(passed),detail=detail))
 print('CHECK',name,bool(passed),flush=True)
def table(name,rows):
 t=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
 t.to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
 return t
def z(a):
 a=np.asarray(a,dtype=float); sd=a.std()
 return (a-a.mean())/sd if sd>1e-14 else a*0
def read10x(path):
 with h5py.File(path) as h:
  g=h['matrix']; X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr().astype(float)
  names=np.array([x.decode() for x in g['features/name'][:]])
  ids=np.array([x.decode() for x in g['features/id'][:]])
  bc=np.array([x.decode() for x in g['barcodes'][:]])
  types=np.array([x.decode() for x in g['features/feature_type'][:]]) if 'feature_type' in g['features'] else np.full(len(names),'unknown')
 return X,names,ids,bc,types
def collapsed(X,names):
 unique,inv=np.unique(names,return_inverse=True)
 C=sparse.csr_matrix((np.ones(len(names)),(np.arange(len(names)),inv)),shape=(len(names),len(unique)))
 return (X@C).tocsr(),unique
def normlog(X):
 umi=np.asarray(X.sum(1)).ravel(); N=(sparse.diags(10000/np.maximum(umi,1))@X).tocsr(); N.data=np.log1p(N.data)
 return N,umi
def graph(xy,radius):
 # Independent dense pairwise distance implementation; small Visium matrices only.
 dist=distance.cdist(xy,xy); A=(dist<=radius+1e-8); np.fill_diagonal(A,False)
 G=sparse.csr_matrix(A,dtype=float); deg=np.asarray(G.sum(1)).ravel()
 return sparse.diags(1/np.maximum(deg,1))@G,deg,G
def effect(y,x,cov,keep):
 v=np.column_stack([y,x]+list(cov)); keep=keep&np.isfinite(v).all(1)
 v=v[keep]; Z=np.column_stack([np.ones(len(v))]+[z(v[:,i]) for i in range(1,v.shape[1])])
 return dict(n=int(keep.sum()),beta=float(np.linalg.lstsq(Z,z(v[:,0]),rcond=None)[0][1]),rank=int(np.linalg.matrix_rank(Z)),parameters=Z.shape[1],condition=float(np.linalg.cond(Z)))

def audit_visium():
 qual=[];geometry=[];models=[];reproductions=[];marker=[];axes=[];coverage=[]
 oldeffects=pd.read_csv(P/'results/mixed_spot_effects_descriptive.tsv',sep='\t')
 folder=B/'expanded_screen_20260925/spatial_data/GSE217414'
 for f in sorted(folder.glob('*_filtered_feature_bc_matrix.h5')):
  sample=f.name.split('_')[0]; sha(f)
  raw,names,ids,bc,types=read10x(f)
  chk(sample+'_count_valid',np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.round(raw.data)).all(),{'feature_types':dict(zip(*np.unique(types,return_counts=True)))})
  X,names=collapsed(raw,names)
  sfolder=next((folder/(sample+'_unpacked')).glob('*spatial'))
  pos=pd.read_csv(sfolder/'tissue_positions_list.csv',header=None,names=['barcode','tissue','ar','ac','py','px']).set_index('barcode')
  sha(sfolder/'tissue_positions_list.csv'); sha(sfolder/'scalefactors_json.json')
  chk(sample+'_barcodes',len(set(bc))==len(bc) and pos.index.is_unique and set(bc)<=set(pos.index),{})
  pos=pos.loc[bc]; sf=json.loads((sfolder/'scalefactors_json.json').read_text())
  umi=np.asarray(X.sum(1)).ravel(); keep=(umi>=500)&pos.tissue.eq(1).to_numpy()
  X=X[keep]; pos=pos.iloc[np.flatnonzero(keep)]; bc=bc[keep]
  N,umi=normlog(X); target=N[:,np.flatnonzero(names=='CGNL1')[0]].toarray().ravel(); targetraw=X[:,np.flatnonzero(names=='CGNL1')[0]].toarray().ravel()
  scores={}
  for key,genes in panels.items():
   ix=np.flatnonzero(np.isin(names,genes)); scores[key]=np.asarray(N[:,ix].mean(1)).ravel()
   coverage.append(dict(sample=sample,panel=key,covered=len(ix),requested=len(genes)))
  saved=pd.read_csv(P/'results'/f'{sample}_spot_profiles.tsv.gz',sep='\t').set_index('barcode').loc[bc]
  maxdiff=max(float(np.max(np.abs(saved[k].to_numpy()-v))) for k,v in scores.items())
  chk(sample+'_profile_reproduced',maxdiff<2e-6 and np.allclose(saved.target_raw,targetraw),{'max_score_abs_difference':maxdiff})
  pix=pos[['px','py']].to_numpy(float)
  ideal=np.c_[pos.ac.to_numpy()*50,pos.ar.to_numpy()*np.sqrt(3)*50]
  xy55=pix*55/sf['spot_diameter_fullres']; xy65=pix*65/sf['spot_diameter_fullres']
  aff=np.linalg.lstsq(np.c_[ideal,np.ones(len(ideal))],pix,rcond=None)[0]
  pixfit=np.c_[ideal,np.ones(len(ideal))]@aff
  oldpitch=float(np.median(cKDTree(xy55).query(xy55,k=2)[0][:,1])); newpitch=float(np.median(cKDTree(xy65).query(xy65,k=2)[0][:,1]))
  geometry.append(dict(sample=sample,spot_diameter_fullres=sf['spot_diameter_fullres'],old_pitch_um=oldpitch,corrected_pitch_um=newpitch,ideal_pitch_um=100,array_fit_rms_pixels=float(np.sqrt(np.mean((pix-pixfit)**2))),old_radius200_actual_um=200*65/55))
  qual.append(dict(sample=sample,n_raw=len(keep),n_QC=len(bc),CGNL1_counts=targetraw.sum(),positive_spots=int((targetraw>0).sum()),positive_fraction=float((targetraw>0).mean()),n_genes=len(names),duplicate_symbols_collapsed=int(raw.shape[1]-len(names))))
  basecov=[scores['endothelial_identity'],np.log1p(umi),scores['fibroblast'],scores['epithelial'],scores['hepatocyte']]
  for geom,xy in [('original55',xy55),('corrected65',xy65),('ideal_lattice100',ideal)]:
   for radius in [100,200,300]:
    W,deg,G=graph(xy,radius); mural=np.asarray(W@scores['mural']).ravel(); tn=np.asarray(W@scores['tcell']).ravel()
    for policy,mask in [('legacy_zero_neighbor',np.ones(len(bc),bool)),('valid_neighbors_only',deg>0)]:
     rr=effect(target,mural,basecov,mask)
     models.append(dict(sample=sample,geometry=geom,radius_um=radius,policy=policy,isolated=int((deg==0).sum()),median_neighbors=float(np.median(deg)),**rr))
     if geom=='original55' and policy=='legacy_zero_neighbor':
      old=float(oldeffects[(oldeffects['sample']==sample)&(oldeffects.radius_um==radius)].beta_neighbor_mural_adjusted.iloc[0])
      reproductions.append(dict(sample=sample,radius=radius,original=old,recomputed=rr['beta'],difference=rr['beta']-old,passed=abs(rr['beta']-old)<2e-6))
    if geom=='corrected65' and radius==200:
     valid=deg>0
     for variant,extra in [('frozen_adjustment',[]),('plus_focal_mural',[scores['mural']])]:
      marker.append(dict(sample=sample,variant=variant,**effect(target,mural,basecov+extra,valid)))
     for excluded in panels['mural']:
      genes=[g for g in panels['mural'] if g!=excluded]; ix=np.flatnonzero(np.isin(names,genes)); v=np.asarray(N[:,ix].mean(1)).ravel()
      marker.append(dict(sample=sample,variant='leave_out_'+excluded,**effect(target,np.asarray(W@v).ravel(),basecov,valid)))
     for name,outcome in [('adhesion',scores['adhesion']),('notch_response',scores['notch_response']),('neighbor_Tcell',tn)]:
      axes.append(dict(sample=sample,endpoint=name,geometry='corrected65',radius_um=200,**effect(outcome,target,basecov,valid),P=np.nan,FDR=np.nan,evidence='descriptive_mixed_spot_no_patient_inference'))
  print('VISIUM_DONE',sample,flush=True)
 table('visium_input_checks',qual);table('visium_geometry',geometry);table('visium_effects_all_fixed_variants',models);table('visium_original_reproduction',reproductions);table('visium_marker_and_local_composition_diagnostics',marker);table('visium_functional_axes_descriptive',axes);table('visium_panel_coverage',coverage)
 chk('all_12_original_visium_coefficients_reproduced',all(r['passed'] for r in reproductions),{})

def audit_ec():
 reg=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text())
 libraries=[];old_profiles=[]
 spec=json.loads((B/'review_consensus_retest_20261002_v1/SPECIFICATION.json').read_text())
 for rec in reg:
  sha(rec['path']);sha(rec['genes']);sha(rec['cells'])
  names=pd.read_csv(rec['genes'],sep='\t').iloc[:,0].astype(str).to_numpy();X=sparse.load_npz(rec['path']).tocsr()
  if X.shape[0]==len(names): X=X.T.tocsr()
  assert X.shape[1]==len(names)
  cellmeta=pd.read_csv(rec['cells'],sep='\t'); assert len(cellmeta)==X.shape[0]
  umi=np.asarray(X.sum(1)).ravel();ng=X.getnnz(1)
  chk(rec['sample']+'_cache_cell_UMI',len(umi)==0 or np.allclose(umi,cellmeta.UMI_QC),{'cells':len(umi),'n_features':len(names)})
  X,names=collapsed(X,names)
  total=np.asarray(X.sum(0)).ravel(); profile=np.log1p(total/max(total.sum(),1)*10000)
  ti=np.flatnonzero(names=='CGNL1'); assert len(ti)==1
  tv=X[:,ti[0]].toarray().ravel()
  rr=dict(dataset=rec['dataset'],sample=rec['sample'],patient=rec['donor'],cells=X.shape[0],total_UMI=total.sum(),CGNL1_count=float(tv.sum()),positive_cells=int((tv>0).sum()),mean_logUMI=float(np.log1p(umi).mean()) if len(umi) else np.nan,empty_sample=X.shape[0]==0)
  for key,genes in {'CGNL1':['CGNL1'],**panels}.items():
   ix=np.flatnonzero(np.isin(names,genes));rr[key]=float(profile[ix].mean()) if len(ix)>=(1 if key=='CGNL1' else 2) and X.shape[0] else np.nan
  libraries.append(rr)
  mt=pd.to_numeric(cellmeta.mito_retest,errors='raise').to_numpy(float)
  for stage in spec['continuous_stages']:
   mask=(umi>=stage['umi'])&(ng>=stage['genes'])&(np.isnan(mt)|(mt<=stage['mt']))
   if stage['branch']=='exclude_reference_immune':mask &= ~cellmeta.reference_agreement.eq('both_immune').to_numpy()
   if stage['branch']=='both_references_endo':mask &= cellmeta.reference_agreement.eq('both_endothelial').to_numpy()
   if mask.sum()<stage['mincells']:continue
   tot=np.asarray(X[mask].sum(0)).ravel();v=np.log1p(tot/max(tot.sum(),1)*1e4)
   old_profiles.append(dict(stage=stage['stage'],dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],CGNL1=v[ti[0]],junction=float(v[np.flatnonzero(np.isin(names,spec['panel']))].mean()),mean_logUMI=float(np.log1p(umi[mask]).mean())))
 lib=table('EC_reconstructed_libraries',libraries)
 oldlib=pd.read_csv(P/'results/accepted_EC_library_measurements.tsv',sep='\t').set_index('sample')
 nonempty=lib[lib.cells>0].set_index('sample');keys=['CGNL1','notch_response','adhesion','hypoxia','historical_junction']
 diff=float(np.max(np.abs(nonempty[keys].to_numpy()-oldlib.loc[nonempty.index,keys].to_numpy())))
 chk('EC_library_panels_reproduced_nonempty',diff<1e-8,{'max_abs_diff':diff,'empty_libraries':lib.loc[lib.empty_sample,'sample'].tolist(),'old_zero_value_in_empty_library_is_not_expression_measurement':True})
 pt=lib[lib.cells>=20].groupby(['dataset','patient'],as_index=False).mean(numeric_only=True)
 count=lib[lib.cells>=20].groupby(['dataset','patient'],as_index=False)[['cells','CGNL1_count','positive_cells']].sum().rename(columns={'cells':'total_cells','CGNL1_count':'total_CGNL1_count','positive_cells':'total_positive_cells'})
 pt=pt.merge(count,on=['dataset','patient']); table('EC_n20_patient_measurements',pt)
 saved=pd.read_csv(P/'results/accepted_EC_patient_measurements.tsv',sep='\t').set_index(['dataset','patient']);got=pt.set_index(['dataset','patient'])
 chk('EC_n20_patient_panels_reproduced',np.allclose(saved[keys],got.loc[saved.index,keys]),{'n':len(pt)})
 oldpt=pd.DataFrame(old_profiles).groupby(['stage','dataset','donor'],as_index=False).mean(numeric_only=True)
 table('A05_raw_reconstructed_patient_profiles',oldpt)
 oldmodels=pd.read_csv(B/'review_consensus_retest_20261002_v1/results/continuous_models.tsv',sep='\t');reps=[]
 for _,r in oldmodels.iterrows():
  if r.status!='estimated': continue
  v=oldpt[oldpt.stage==r.stage].copy()
  if r.scope!='pooled': v=v[v.dataset==r.scope]
  F=pd.DataFrame({'CGNL1':v.CGNL1})
  if r.model!='unadjusted':F=pd.concat([F,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
  if r.model=='study_depth_adjusted':F['depth']=v.mean_logUMI
  fit=sm.OLS(v.junction,sm.add_constant(F,has_constant='add')).fit(use_t=True)
  beta=float(fit.params.CGNL1)*v.CGNL1.std()/v.junction.std();pv=float(fit.pvalues.CGNL1)
  reps.append(dict(stage=r.stage,scope=r.scope,model=r.model,family=r.family,n=len(v),beta=beta,P=pv,old_beta=r.beta,old_P=r.P,old_FDR=r.FDR,pass_beta=abs(beta-r.beta)<1e-8,pass_P=abs(pv-r.P)<1e-8))
 rep=pd.DataFrame(reps);rep['FDR']=np.nan
 for fam,ix in rep.groupby('family').groups.items():rep.loc[ix,'FDR']=multipletests(rep.loc[ix,'P'],method='fdr_bh')[1]
 rep['pass_FDR']=np.abs(rep.FDR-rep.old_FDR)<1e-8;table('A05_full_model_reproduction',rep)
 chk('A05_all_estimable_coefficients_P_FDR_reproduced',rep[['pass_beta','pass_P','pass_FDR']].to_numpy().all(),{'models':len(rep)})
 # New exploratory audit on pre-existing exposed n20 patients, not external validation.
 outcomes=['notch_response','adhesion','historical_junction'];models=[];pred=[];loo=[]
 for scope in ['pooled','GSE245552','GSE285991']:
  v=pt if scope=='pooled' else pt[pt.dataset==scope];v=v.reset_index(drop=True)
  base=pd.DataFrame({'intercept':np.ones(len(v)),'mean_logUMI':v.mean_logUMI})
  base=pd.concat([base,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
  F=base.copy();F['CGNL1']=v.CGNL1
  for endpoint in outcomes:
   y=v[endpoint]; non=sm.OLS(y,F).fit();rob=sm.OLS(y,F).fit(cov_type='HC3',use_t=True)
   for kind,fit in [('OLS',non),('HC3',rob)]:
    models.append(dict(scope=scope,endpoint=endpoint,variance=kind,family='AUDIT_EC_'+scope+'_'+kind+'_three_axes',n=len(v),beta=float(fit.params.CGNL1)*v.CGNL1.std()/y.std(),slope=float(fit.params.CGNL1),CI_low=float(fit.conf_int().loc['CGNL1',0]),CI_high=float(fit.conf_int().loc['CGNL1',1]),P=float(fit.pvalues.CGNL1),df=float(fit.df_resid),interpretation='exploratory_existing_exposed_patients_partial_H3_only'))
   for i in range(len(v)):
    mask=np.arange(len(v))!=i
    bf=np.linalg.lstsq(base.to_numpy()[mask],y.to_numpy()[mask],rcond=None)[0];ff=np.linalg.lstsq(F.to_numpy()[mask],y.to_numpy()[mask],rcond=None)[0]
    pred.append(dict(scope=scope,endpoint=endpoint,patient=v.patient.iloc[i],observed=float(y.iloc[i]),baseline=float(base.to_numpy()[i]@bf),with_CGNL1=float(F.to_numpy()[i]@ff)))
    sub=sm.OLS(y[mask],F.loc[mask]).fit()
    loo.append(dict(scope=scope,endpoint=endpoint,omitted=v.patient.iloc[i],slope=float(sub.params.CGNL1),P=float(sub.pvalues.CGNL1),purpose='influence_diagnostic_not_new_discovery_family'))
 mt=pd.DataFrame(models);mt['FDR']=np.nan
 for fam,ix in mt.groupby('family').groups.items():mt.loc[ix,'FDR']=multipletests(mt.loc[ix,'P'],method='fdr_bh')[1]
 table('EC_exploratory_axes_models',mt);pr=table('EC_leave_one_patient_predictions',pred);table('EC_leave_one_patient_influence',loo)
 perf=[]
 for (scope,endpoint),g in pr.groupby(['scope','endpoint']):
  bse=(g.observed-g.baseline)**2;ase=(g.observed-g.with_CGNL1)**2
  perf.append(dict(scope=scope,endpoint=endpoint,n=len(g),MSE_baseline=bse.mean(),MSE_with_CGNL1=ase.mean(),relative_MSE_improvement=1-ase.mean()/bse.mean(),patients_improved=int((ase<bse).sum()),baseline='study+mean_logUMI',missing_full_H3_covariates='endothelial abundance/subtype/stress; existing patients, no independent cohort'))
 table('EC_incremental_prediction_diagnostic',perf)

def audit_m242():
 ap=B/'annotation_panel_20261007_v1'; fp=B/'fullslice_20261007_v1'
 rawpath=B/'spatial_v3/raw/GSM9815454_tumor1_M242T_cell_feature_matrix.h5';sha(rawpath)
 X,names,ids,bc,types=read10x(rawpath)
 aliases={'CGNL1','JACOP','PCING','ENSG00000128849','84952'}
 chk('M242T_target_absence_confirmed',not aliases.intersection(set(names)|set(ids)),{'features':len(names),'gene_expression_features':int((types=='Gene Expression').sum())})
 a=ad.read_h5ad('/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1/GSM9815454_panel_annotated_v1.h5ad',backed='r')
 obs=a.obs.copy();xy=np.array(a.obsm['spatial']);a.file.close()
 with gzip.open(BIG/'downloads/GSM9815454_cells.parquet.gz','rb') as f:original=pd.read_parquet(io.BytesIO(f.read())).set_index('cell_id')
 original=original.loc[bc];obs=obs.loc[bc]
 chk('M242T_raw_counts_metadata_alignment',np.array_equal(np.asarray(X.sum(1)).ravel(),original.transcript_counts.to_numpy()) and np.array_equal(original.index,obs.index),{'n':len(bc)})
 chk('M242T_coordinates_alignment',np.allclose(xy,original[['x_centroid','y_centroid']]),{})
 q=(original.transcript_counts>=50)&(original.nucleus_count==1)
 chk('M242T_QC_reproduction',np.array_equal(q.to_numpy(),obs.qc_pass.to_numpy()),{'passed':int(q.sum())})
 obs=obs.loc[q];xy=xy[q.to_numpy()];labels=obs.annotation.astype(str).to_numpy();classes=sorted(set(labels));codes=pd.Categorical(labels,categories=classes).codes
 labelhash=hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(obs.index,labels)).encode()).hexdigest();expected=json.loads((fp/'standards/FULL_LABEL_FREEZE.json').read_text())['labels_sha256']
 chk('M242T_labels_frozen',labelhash==expected,{'sha256':labelhash})
 choose=np.random.default_rng(20261008).choice(len(obs),256,replace=False);indices=list(obs.index[choose]);rows=[]
 for cond in ['knn_15','knn_30','knn_60','radius_100']:
  old=pd.read_csv(fp/'results'/f'{cond}_composition.tsv.gz',sep='\t',index_col=0).loc[indices,classes]
  maximum=0.;mismatch=0
  for j,i in enumerate(choose):
   d=np.linalg.norm(xy-xy[i],axis=1);d[i]=np.inf
   nn=np.argsort(d)[:int(cond.split('_')[1])] if cond.startswith('knn') else np.flatnonzero(d<100)
   result=np.bincount(codes[nn],minlength=len(classes))/len(nn) if len(nn)>=3 else np.full(len(classes),np.nan)
   if not np.allclose(result,old.iloc[j],equal_nan=True):mismatch+=1
   if np.isfinite(result).all():maximum=max(maximum,float(np.abs(result-old.iloc[j]).max()))
  rows.append(dict(condition=cond,audited_centers=len(choose),mismatch=mismatch,max_abs_difference=maximum))
 table('M242T_independent_neighborhood_reproduction',rows);chk('M242T_1024_profiles_reproduced',sum(x['mismatch'] for x in rows)==0,{})
 reviewed=pd.read_csv(fp/'results/fullslice_cell_neighborhoods.tsv.gz',sep='\t',index_col=0)
 ari=adjusted_rand_score(reviewed.primary_niche,reviewed.alternate_seed_cluster)
 chk('M242T_discrete_seed_ARI_reproduced',abs(ari-json.loads((fp/'results/ANALYSIS_ACCEPTANCE.json').read_text())['seed_ARI'])<1e-10,{'ARI':ari})
 table('M242T_annotation_counts_audited',pd.Series(labels).value_counts().rename_axis('label').reset_index(name='cells'))
 # Read-only held-out vascular marker support; not reference-derived accuracy.
 XX=X[q.to_numpy()];gene={g:i for i,g in enumerate(names)};support=[]
 for label in ['endothelial','mural','fibroblast','hepatocyte_lineage','uncertain']:
  mask=labels==label
  for gn in ['PECAM1','VWF','KDR','FLT1','ROBO4','RGS5','ACTA2','TAGLN','EPCAM','KRT19','CYP2E1']:
   if gn in gene:
    v=XX[mask,gene[gn]].toarray().ravel();support.append(dict(label=label,gene=gn,cells=len(v),detection_fraction=float(np.mean(v>0)),mean_raw=float(v.mean())))
 table('M242T_marker_support_audited',support)

def audit_methods():
 inputs=json.loads((P/'results/TECHNICAL_INPUTS.json').read_text());samplegraphs=[]
 for rec in inputs:
  a=ad.read_h5ad(rec['path']); G=a.obsp['spatial_connectivities']
  samplegraphs.append(G)
  expected=pd.read_csv(P/'results'/f'{rec["sample"]}_spot_profiles.tsv.gz',sep='\t').set_index('barcode')
  chk(rec['sample']+'_technical_object',a.obs_names.is_unique and set(a.obs.barcode)<=set(expected.index) and G.diagonal().sum()==0,{'spots':a.n_obs,'edges':G.nnz})
  sc=ad.read_h5ad(BIG/'inputs'/f'{rec["sample"]}_LIANA_scores.h5ad')
  dense=sc.X.toarray() if sparse.issparse(sc.X) else np.asarray(sc.X)
  tab=pd.read_csv(P/'results'/f'{rec["sample"]}_LIANA_pair_scores.tsv',sep='\t').set_index('pair')
  chk(rec['sample']+'_LIANA_saved_score_means',np.isfinite(dense).all() and np.allclose(dense.mean(0),tab.loc[sc.var_names,'mean_score']),{'pairs':sc.n_vars,'cell_specific_sender_receiver':False,'input_geometry_inherits_55_conversion':True})
 cfg=json.loads((P/'results/NICHECOMPASS_PARAMETERS.json').read_text());required={'CGNL1'}|set(sum([panels[k] for k in ['notch_response','adhesion','historical_junction']],[]))
 chk('NicheCompass_exclusion_specification',required<=set(cfg['target_exclusions']),{'required':sorted(required)})
 a=pd.read_csv(P/'results/NicheCompass_clusters_20261006.tsv',sep='\t');b=pd.read_csv(P/'results/NicheCompass_clusters_20261007.tsv',sep='\t').set_index('spot').loc[a.spot]
 ari=adjusted_rand_score(a.technical_cluster,b.technical_cluster)
 chk('NicheCompass_saved_ARI_reproduced',abs(ari-json.loads((P/'results/NICHECOMPASS_COMPLETE.json').read_text())['cross_seed_cluster_ARI'])<1e-10,{'ARI':ari})
 for seed in [20261006,20261007]:
  f=P/'results'/f'NicheCompass_latent_{seed}.tsv.gz'; lat=pd.read_csv(f,sep='\t',index_col=0)
  chk('NicheCompass_latent_'+str(seed),len(lat)==3000 and np.isfinite(lat.to_numpy()).all(),{'shape':list(lat.shape)})
 # Reconstruct actual chosen feature list, independently from masks and saved parameters.
 arr=[ad.read_h5ad(r['path']) for r in inputs];common=arr[0].var_names.intersection(arr[1].var_names)
 X=sparse.vstack([a[:,common].X for a in arr]).tocsr();N,u=normlog(X);vv=np.asarray(N.power(2).mean(0)).ravel()-np.asarray(N.mean(0)).ravel()**2
 chosen=[str(common[i]) for i in np.argsort(-vv,kind='stable') if str(common[i]) not in required][:800]
 selected=set(chosen)|set(g for pair in cfg['prior_programs'] for g in pair.split('_'))
 chk('NicheCompass_reconstructed_input_excludes_targets',not required.intersection(selected),{'features':len(selected),'declared_features':cfg['genes']})

for name,fun in [('visium',audit_visium),('endothelial',audit_ec),('M242T',audit_m242),('methods',audit_methods)]:
 print('MODULE_START',name,flush=True)
 try:
  fun();modules[name]={'execution':'complete'}
 except Exception as exc:
  modules[name]={'execution':'failed','error':repr(exc),'traceback':traceback.format_exc()};print(traceback.format_exc(),flush=True)
 (OUT/'MODULE_STATUS.json').write_text(json.dumps(modules,indent=2))
 (OUT/'CHECKS.json').write_text(json.dumps(checks,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
 (OUT/'SOURCE_FINGERPRINTS.json').write_text(json.dumps(fingerprints,indent=2))
done={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'modules':modules,'checks_passed':sum(x['passed'] for x in checks),'checks_failed':sum(not x['passed'] for x in checks),'original_analysis_files_modified':False,'audit_scope_sha256':sha(R/'AUDIT_SCOPE.json')}
(OUT/'NUMERIC_AUDIT_COMPLETE.json').write_text(json.dumps(done,indent=2));print(json.dumps(done),flush=True)
raise SystemExit(0 if all(v['execution']=='complete' for v in modules.values()) else 1)
