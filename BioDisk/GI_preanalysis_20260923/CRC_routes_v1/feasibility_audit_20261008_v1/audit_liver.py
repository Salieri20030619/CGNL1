from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,tarfile,gzip,io,datetime,hashlib
import numpy as np,pandas as pd
from scipy import sparse,io as sio,stats
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'feasibility_audit_20261008_v1';O=R/'results'
P=json.loads((B/'v3_pilot_20261006/plan_v3.json').read_text());panels=P['panels'];needed=sorted(set(['CGNL1']+sum(panels.values(),[])))
jobs=json.loads((B/'spatial_v3/metadata/LIVER_DOWNLOAD_PLAN.json').read_text());ann=pd.read_csv(B/'results/spatial_author_metadata.tsv.gz',sep='\t');ann['sample']=ann.Barcode.str.split('~').str[0];ann['barcode']=ann.Barcode.str.split('~').str[1]
legacy=pd.read_csv(B/'spatial_v3/results/screen/spot_target_coordinates.tsv.gz',sep='\t');checks=[];models=[];covrows=[];integrity=[]
def z(x):
 x=np.asarray(x,float);return (x-x.mean())/x.std() if x.std()>1e-12 else x*0
def fit(y,x,cov,mask,region=None):
 Z=np.column_stack([x]+list(cov));valid=mask&np.isfinite(Z).all(1)&np.isfinite(y)
 if valid.sum()<20:return dict(n=int(valid.sum()),beta=np.nan,status='insufficient_spots')
 F=np.column_stack([np.ones(valid.sum())]+[z(v) for v in Z[valid].T])
 if region is not None:F=np.c_[F,pd.get_dummies(pd.Series(region[valid]),drop_first=True,dtype=float).to_numpy()]
 rank=np.linalg.matrix_rank(F);cond=np.linalg.cond(F)
 if rank<F.shape[1] or cond>1e8:return dict(n=int(valid.sum()),beta=np.nan,status='rank_or_condition',rank=int(rank),condition=float(cond))
 return dict(n=int(valid.sum()),beta=float(np.linalg.lstsq(F,z(y[valid]),rcond=None)[0][1]),status='estimated_descriptive',rank=int(rank),condition=float(cond))
for job in jobs:
 name=job['sample'];path=B/'spatial_v3/raw/liver'/f'{name}.tar.gz';print('LIVER_START',name,flush=True)
 files={};X=None
 # Streaming traversal avoids extracting large images or repeatedly decompressing the same archive.
 with tarfile.open(path,'r|gz') as tar:
  for m in tar:
   if not m.isfile():continue
   if m.name.endswith('filtered_feature_bc_matrix/features.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as f:features=pd.read_csv(f,sep='\t',header=None)
   elif m.name.endswith('filtered_feature_bc_matrix/barcodes.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as f:barcodes=pd.read_csv(f,sep='\t',header=None)[0].astype(str).to_numpy()
   elif m.name.endswith('filtered_feature_bc_matrix/matrix.mtx.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(m)) as f,threadpool_limits(limits=4):X=sio.mmread(f).tocsr()
   elif m.name.endswith('spatial/tissue_positions.csv'):files['pos']=tar.extractfile(m).read()
   elif m.name.endswith('spatial/tissue_positions_list.csv'):files['poslist']=tar.extractfile(m).read()
 names=features[1].astype(str).to_numpy();assert X.shape==(len(names),len(barcodes));X.eliminate_zeros()
 a=ann[ann['sample']==name].copy().reset_index(drop=True);ii=pd.Index(barcodes).get_indexer(a.barcode);assert (ii>=0).all();Y=X[:,ii];del X
 umi=np.asarray(Y.sum(0)).ravel();ng=Y.getnnz(0)
 positions=pd.read_csv(io.BytesIO(files['pos']),header=0) if 'pos' in files else pd.read_csv(io.BytesIO(files['poslist']),header=None,names=['barcode','in_tissue','array_row','array_col','py','px'])
 positions=positions.set_index('barcode').loc[a.barcode]
 xy=np.c_[positions.array_col.to_numpy()*50,positions.array_row.to_numpy()*50*np.sqrt(3)]
 saved=legacy[legacy['sample']==name].set_index('barcode').loc[a.barcode]
 target=np.asarray(Y[np.flatnonzero(names=='CGNL1')].sum(0)).ravel()
 good=(umi>=500)&(ng>=200)
 rec='Recurrence' in job['sample_type'];base=dict(sample=name,patient=job['patient'],recurrence=rec)
 match=np.array_equal(umi,saved.UMI)&np.array_equal(target,saved.CGNL1_counts)&np.array_equal(good,saved.QC_pass)&np.allclose(xy,saved[['x_um','y_um']])
 integrity.append({**base,'author_spots':len(a),'QC_spots':int(good.sum()),'target_counts_QC':float(target[good].sum()),'raw_vs_saved_pass':bool(match),'target_feature_rows':int((names=='CGNL1').sum())})
 assert match,(name,'source mismatch')
 # Direct distance reconstruction from author pathology labels for the historical screen.
 hepatic=a.Layer3.eq('Hepatic lobule').to_numpy();distance=cKDTree(xy[hepatic]).query(xy)[0] if hepatic.any() else np.full(len(xy),np.nan)
 checks.append(dict(sample=name,distance_reproduced=bool(np.allclose(distance,saved.distance_to_hepatic_um,equal_nan=True)),annotation_matched=bool(np.array_equal(a.Layer3,saved.Layer3))))
 # Scores computed from measured genes after symbol summation, preserving panel availability.
 values={}
 for gn in needed:
  ix=np.flatnonzero(names==gn)
  if len(ix):values[gn]=np.log1p(np.asarray(Y[ix].sum(0)).ravel()/np.maximum(umi,1)*10000)
 scores={}
 for key,gs in panels.items():
  available=[g for g in gs if g in values];covrows.append({**base,'panel':key,'covered':len(available),'expected':len(gs),'missing':'|'.join(g for g in gs if g not in values)})
  scores[key]=np.mean([values[g] for g in available],axis=0) if len(available)>=2 else np.full(len(a),np.nan)
 del Y
 idx=np.flatnonzero(good);xq=xy[idx];sq={k:v[idx] for k,v in scores.items()};t=values['CGNL1'][idx];uq=umi[idx];regions=a.Layer3.to_numpy()[idx]
 cov=[sq['endothelial_identity'],np.log1p(uq),sq['fibroblast'],sq['epithelial'],sq['hepatocyte']]
 tree=cKDTree(xq)
 for radius in [100,200,300]:
  pp=tree.query_pairs(radius+1e-6,output_type='ndarray');row=np.r_[pp[:,0],pp[:,1]];col=np.r_[pp[:,1],pp[:,0]];G=sparse.csr_matrix((np.ones(len(row)),(row,col)),shape=(len(idx),len(idx)))
  deg=np.asarray(G.sum(1)).ravel();W=sparse.diags(1/np.maximum(deg,1))@G;nm=np.asarray(W@sq['mural']).ravel();nt=np.asarray(W@sq['tcell']).ravel()
  for branch in ['all_QC','all_QC_plus_author_region_intercepts','macro_tumor_QC_only']:
   mask=deg>0
   if branch=='macro_tumor_QC_only':mask &= regions=='Liver macrometastasis tumor'
   reg=regions if branch=='all_QC_plus_author_region_intercepts' else None
   for hypothesis,endpoint,y,x in [('H1','neighbor_mural',t,nm),('H2','adhesion',sq['adhesion'],t),('H2','neighbor_Tcell',nt,t),('H4_readout','notch_response',sq['notch_response'],t)]:
    models.append({**base,'radius_um':radius,'branch':branch,'hypothesis':hypothesis,'endpoint':endpoint,**fit(y,x,cov,mask,reg)})
 pd.DataFrame(integrity).to_csv(O/'GSE294385_raw_reproduction.tsv',sep='\t',index=False)
 pd.DataFrame(models).to_csv(O/'GSE294385_v3_section_effects.tsv',sep='\t',index=False)
 print('LIVER_DONE',name,len(idx),flush=True)
pd.DataFrame(checks).to_csv(O/'GSE294385_geometry_annotation_checks.tsv',sep='\t',index=False);pd.DataFrame(covrows).to_csv(O/'GSE294385_v3_panel_coverage.tsv',sep='\t',index=False)
m=pd.DataFrame(models);pt=[];summ=[]
for scope in ['baseline','all_timepoints']:
 v=m[~m.recurrence] if scope=='baseline' else m
 for key,g in v.groupby(['patient','radius_um','branch','hypothesis','endpoint']):
  val=g.beta.dropna()
  if len(val):pt.append(dict(zip(['patient','radius_um','branch','hypothesis','endpoint'],key),scope=scope,beta=float(val.mean()),sections=len(val)))
pp=pd.DataFrame(pt);pp.to_csv(O/'GSE294385_v3_patient_effects.tsv',sep='\t',index=False)
for key,g in pp.groupby(['scope','radius_um','branch','hypothesis','endpoint']):
 v=g.beta.to_numpy();n=len(v);d=dict(zip(['scope','radius_um','branch','hypothesis','endpoint'],key));d.update(n_patients=n,positive=int((v>0).sum()),negative=int((v<0).sum()),mean_beta=float(v.mean()),median_beta=float(np.median(v)),P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=4 and v.std()>0:
  d['P']=float(stats.ttest_1samp(v,0).pvalue);ci=stats.t.interval(.95,n-1,loc=v.mean(),scale=stats.sem(v));d.update(CI_low=float(ci[0]),CI_high=float(ci[1]))
 d['family']='AUDIT_GSE294385_'+d['scope']+'_'+d['hypothesis'];d['evidence']='exploratory_existing_cohort_mixed_spots';summ.append(d)
s=pd.DataFrame(summ);s['FDR']=np.nan
for fam,ix in s.groupby('family').groups.items():
 valid=s.loc[ix].index[s.loc[ix,'P'].notna()]
 if len(valid):s.loc[valid,'FDR']=multipletests(s.loc[valid,'P'],method='fdr_bh')[1]
s.to_csv(O/'GSE294385_v3_patient_summary.tsv',sep='\t',index=False)
# Historical effects reconstructed from validated per-spot raw counts, independent of stored section coefficients.
oldsec=pd.read_csv(B/'spatial_v3/results/screen/section_effects.tsv',sep='\t');reps=[]
for _,r in oldsec[oldsec.gene.eq('CGNL1')].iterrows():
 a=legacy[legacy['sample']==r['sample']].copy();qc=a.QC_pass.to_numpy();mac=a.Layer3.eq('Liver macrometastasis tumor').to_numpy();dep=a.distance_to_hepatic_um.to_numpy()
 if r.contrast.startswith('tumor_'):
  import re
  near,far=map(int,re.findall(r'\d+',r.contrast));m1=qc&mac&(dep<=near+1e-6);m0=qc&mac&(dep>=far-1e-6)
 else:
  label={'hepatic_lobule_vs_macro_tumor':'Hepatic lobule','macro_stroma_vs_macro_tumor':'Liver macrometastasis stroma','micro_tumor_vs_macro_tumor':'Liver micrometastasis tumor'}[r.contrast];m1=qc&a.Layer3.eq(label).to_numpy();m0=qc&mac
 if r.method=='pseudobulk':v=np.log2(1+a.loc[m1,'CGNL1_counts'].sum()/a.loc[m1,'UMI'].sum()*1e6)-np.log2(1+a.loc[m0,'CGNL1_counts'].sum()/a.loc[m0,'UMI'].sum()*1e6)
 elif r.method=='spot_mean':v=a.loc[m1,'CGNL1_log2CPM'].mean()-a.loc[m0,'CGNL1_log2CPM'].mean()
 else:
  use=m1|m0;cov=[a.loc[use,k+'_score'].to_numpy() for k in ['hepatic','stromal','epithelial','endothelial']]+[np.log1p(a.loc[use,'UMI'].to_numpy())];F=np.c_[np.ones(use.sum()),m1[use].astype(float),np.column_stack([z(x) for x in cov])];v=np.linalg.lstsq(F,a.loc[use,'CGNL1_log2CPM'],rcond=None)[0][1]
 reps.append(dict(sample=r['sample'],contrast=r.contrast,method=r.method,recomputed=float(v),old_effect=r.effect,passed=bool(np.isclose(v,r.effect,equal_nan=True))))
pd.DataFrame(reps).to_csv(O/'GSE294385_historical_CGNL1_section_reproduction.tsv',sep='\t',index=False)
oldpt=pd.read_csv(B/'spatial_v3/results/screen/patient_effects.tsv',sep='\t');oldsum=pd.read_csv(B/'spatial_v3/results/screen/summary.tsv',sep='\t');rr=[]
for _,row in oldsum.iterrows():
 g=oldpt[(oldpt.scope==row.scope)&(oldpt.gene==row.gene)&(oldpt.contrast==row.contrast)&(oldpt.method==row.method)];v=g.effect.to_numpy();pv=float(stats.ttest_1samp(v,0).pvalue) if len(v)>=10 else np.nan
 rr.append(dict(scope=row.scope,gene=row.gene,contrast=row.contrast,method=row.method,n=len(v),mean_effect=float(v.mean()),P=pv,old_P=row.P,old_FDR=row.FDR,family=row.family,pass_effect=bool(np.isclose(v.mean(),row.mean_effect)),pass_P=bool(np.isclose(pv,row.P,equal_nan=True))))
pd.DataFrame(rr).to_csv(O/'GSE294385_historical_patient_reproduction.tsv',sep='\t',index=False)
done={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archives':len(integrity),'mapped_patients':len(set(x['patient'] for x in integrity)),'raw_all24_match':all(x['raw_vs_saved_pass'] for x in integrity),'geometry_all_match':all(x['distance_reproduced'] and x['annotation_matched'] for x in checks),'historical_CGNL1_section_effects_reproduced':all(x['passed'] for x in reps),'historical_patient_effects_and_P_reproduced':all(x['pass_effect'] and x['pass_P'] for x in rr),'new_analysis':'exploratory_frozen_audit_branches_no_spot_inference'}
(O/'LIVER_AUDIT_COMPLETE.json').write_text(json.dumps(done,indent=2));print(json.dumps(done),flush=True)
