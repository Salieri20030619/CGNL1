from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,hashlib,itertools
import numpy as np,pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests
R=Path.cwd();S=json.loads((R/'ANALYSIS_SPECIFICATION.json').read_text());F=json.loads((R/'SPEC_FREEZE.json').read_text());assert hashlib.sha256((R/'ANALYSIS_SPECIFICATION.json').read_bytes()).hexdigest()==F['sha256']
assert (R/'EXTRACTION_COMPLETE.json').exists()
amend={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'target_effects_not_yet_seen':True,'implementation_details':{'Layer2':'Liver_metastatic carcinoma','block_sensitivity':'250um square tiles within same macro component and near/deep group; mean measurements; >=5 spots per retained tile; >=3 tiles per band; no tile P values','neighbor_adjustment':'add 200um neighbor epithelial, fibroblast, endothelial and hepatocyte scores to main covariates','matching':'greedy global nearest available pair without replacement within section/component, RMS z-distance <=0.75; >=20 pairs/section; report covariate SMD and flag max>0.2','rank':'drop only constant nuisance columns, require full rank and residual df>=10; retain failure rows','inference':'exact sign-flip with two-sided mean statistic; t CI descriptive; sign-test P sensitivity; no inference for n<2; all n<10 labeled exploratory'}}
(R/'IMPLEMENTATION_FREEZE.json').write_text(json.dumps(amend,indent=2),encoding='utf-8')
d=pd.read_csv(R/'results/spatial_spot_features.tsv.gz',sep='\t');targets=S['targets'];covs=S['main_covariates'];rows=[];elig=[];checks=[];support=[];matched=[]
variants=[('primary',50,200,500,1500,200,covs),('depth_only',50,200,500,1500,200,['logUMI']),('bands_alt',100,300,600,1500,200,covs),('radius100',50,200,500,1500,100,covs),('radius300',50,200,500,1500,300,covs),('neighbor_adjusted',50,200,500,1500,200,covs+['neighbor_epithelial_identity','neighbor_fibroblast','neighbor_endothelial','neighbor_hepatocyte']),('blocks250',50,200,500,1500,200,covs)]
def select(v,low,high,deep,maxdeep):
 v=v[v.region.eq('Liver macrometastasis tumor')&v.Layer2.eq('Liver_metastatic carcinoma')].copy()
 v['near']=v.distance.between(low,high).astype(float);far=v.distance.between(deep,maxdeep);v=v[v.near.eq(1)|far].copy()
 good=[]
 for k,g in v.groupby('macro_component'):
  if g.near.sum()>=5 and (1-g.near).sum()>=5:good.append(k)
 return v[v.macro_component.isin(good)].copy()
def design(v,cc,radius,endpoint):
 Z=v[cc].to_numpy(float);mu=Z.mean(0);sd=Z.std(0);active=sd>1e-10
 component=pd.get_dummies(v.macro_component,drop_first=True,dtype=float).to_numpy();base=np.c_[np.ones(len(v)),(Z[:,active]-mu[active])/sd[active],component]
 near=v.near.to_numpy();names=['intercept']+list(np.array(cc)[active])+['component'+str(i) for i in range(component.shape[1])]
 if endpoint=='E1':return np.c_[base,near],len(names),names+['near']
 nm=v['neighbor_mural_'+str(radius)].to_numpy();ns=nm.std()
 if ns<=1e-10:return None,None,None
 nz=(nm-nm.mean())/ns;return np.c_[base,near,nz,near*nz],len(names)+2,names+['near','neighbor_mural','near_x_neighbor_mural']
for sample,v0 in d.groupby('sample',sort=True):
 patient=int(v0.patient.iloc[0]);print('MODEL',sample,flush=True)
 for variant,lo,hi,deep,maxdeep,radius,cc in variants:
  v=select(v0,lo,hi,deep,maxdeep);v=v.replace([np.inf,-np.inf],np.nan).dropna(subset=cc+targets+['neighbor_mural_'+str(radius)])
  nnear=int(v.near.sum());nfar=len(v)-nnear;info={'sample':sample,'patient':patient,'variant':variant,'near_spots':nnear,'deep_spots':nfar}
  if nnear<20 or nfar<20:
   elig.append(dict(info,endpoint='both',status='insufficient_within_component_band_coverage'));continue
  if variant=='primary':
   for band,g in v.groupby('near'):
    support.append(dict(sample=sample,patient=patient,band='near' if band else 'deep',spots=len(g),CGNL1_detection=float((g.CGNL1_count>0).mean()),CGNL1_counts=int(g.CGNL1_count.sum()),mean_expression=float(g.CGNL1.mean()),mean_distance=float(g.distance.mean()),components=int(g.macro_component.nunique())))
   # Target-blind composition matching; outcome never enters distance or ordering.
   Z=v[covs].to_numpy();sd=Z.std(0);active=sd>1e-10;Z=(Z[:,active]-Z[:,active].mean(0))/sd[active];pairs=[]
   for comp in sorted(v.macro_component.unique()):
    a=np.flatnonzero(v.near.eq(1).to_numpy()&v.macro_component.eq(comp).to_numpy());b=np.flatnonzero(v.near.eq(0).to_numpy()&v.macro_component.eq(comp).to_numpy())
    dm=cdist(Z[a],Z[b])/np.sqrt(active.sum());order=np.argsort(dm,axis=None,kind='stable');useda=set();usedb=set()
    for t in order:
     ia,ib=np.unravel_index(t,dm.shape)
     if dm[ia,ib]>.75:break
     if ia not in useda and ib not in usedb:useda.add(ia);usedb.add(ib);pairs.append((a[ia],b[ib],float(dm[ia,ib])))
   if pairs:
    aa=np.array([p[0] for p in pairs]);bb=np.array([p[1] for p in pairs]);smd=np.mean(Z[aa]-Z[bb],axis=0);matched.append(dict(sample=sample,patient=patient,pairs=len(pairs),max_abs_SMD=float(abs(smd).max()),balanced=bool(abs(smd).max()<=.2),mean_distance=float(np.mean([p[2] for p in pairs]))))
    pd.DataFrame([dict(sample=sample,patient=patient,near_barcode=v.iloc[a].barcode,deep_barcode=v.iloc[b].barcode,distance=dis) for a,b,dis in pairs]).to_csv(R/'results'/('matched_pairs_'+sample+'.tsv'),sep='\t',index=False)
    if len(pairs)>=20:
     for gene in targets:
      y=v[gene].to_numpy();effect=float((y[aa]-y[bb]).mean());ys=y.std()
      rows.append(dict(sample=sample,patient=patient,variant='matched',endpoint='E1',gene=gene,effect=effect,effect_SD=effect/ys if ys>1e-10 else np.nan,n=len(pairs),balanced=bool(abs(smd).max()<=.2),condition=np.nan))
  if variant=='blocks250':
   v=v.assign(tilex=np.floor(v.x_um/250).astype(int),tiley=np.floor(v.y_um/250).astype(int));gr=v.groupby(['macro_component','near','tilex','tiley']);cnt=gr.size();v=gr[[*cc,*targets,'neighbor_mural_'+str(radius)]].mean().loc[cnt>=5].reset_index();nnear=int(v.near.sum());nfar=len(v)-nnear
   if min(nnear,nfar)<3:elig.append(dict(info,endpoint='both',status='insufficient_spatial_blocks'));continue
  for endpoint in ['E1','E2']:
   D,index,names=design(v,cc,radius,endpoint)
   if D is None:elig.append(dict(info,endpoint=endpoint,status='zero_neighbor_variance'));continue
   rank=np.linalg.matrix_rank(D);cond=np.linalg.cond(D)
   if rank<D.shape[1] or len(v)-D.shape[1]<10 or cond>1e8:elig.append(dict(info,endpoint=endpoint,status='rank_df_condition_failure',rank=rank,parameters=D.shape[1],n=len(v),condition=cond));continue
   elig.append(dict(info,endpoint=endpoint,status='estimated',rank=rank,parameters=D.shape[1],n=len(v),condition=cond))
   Y=v[targets].to_numpy();beta=np.linalg.lstsq(D,Y,rcond=None)[0][index];N=np.delete(D,index,axis=1);x=D[:,index];rx=x-N@np.linalg.lstsq(N,x,rcond=None)[0];ry=Y-N@np.linalg.lstsq(N,Y,rcond=None)[0];b2=rx@ry/(rx@rx)
   checks.append({'sample':sample,'variant':variant,'endpoint':endpoint,'max_FWL_error':float(np.max(abs(beta-b2)))})
   for gene,b,ys in zip(targets,beta,Y.std(0)):
    rows.append(dict(sample=sample,patient=patient,variant=variant,endpoint=endpoint,gene=gene,effect=float(b),effect_SD=float(b/ys) if ys>1e-10 else np.nan,n=len(v),condition=float(cond)))
q=pd.DataFrame(rows);q.to_csv(R/'results/section_effects.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(R/'results/model_eligibility.tsv',sep='\t',index=False);pd.DataFrame(support).to_csv(R/'results/target_measurement_support.tsv',sep='\t',index=False);pd.DataFrame(matched).to_csv(R/'results/matching_balance.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(R/'results/FWL_checks.tsv',sep='\t',index=False)
p=q.groupby(['variant','endpoint','gene','patient'],as_index=False).agg(effect=('effect','mean'),effect_SD=('effect_SD','mean'),sections=('sample','nunique'));p.to_csv(R/'results/patient_effects.tsv',sep='\t',index=False)
def infer(x):
 x=np.asarray(x,float);x=x[np.isfinite(x)];n=len(x);out={'n':n,'mean':float(x.mean()) if n else np.nan,'median':float(np.median(x)) if n else np.nan,'positive':int((x>0).sum()),'negative':int((x<0).sum()),'P':np.nan,'sign_P':np.nan,'CI_low':np.nan,'CI_high':np.nan,'LOO_mean_min':np.nan,'LOO_mean_max':np.nan,'evidence':'exploratory_small_n' if n<10 else 'exploratory_patient_level'}
 if n>=2:
  signs=np.array(list(itertools.product([-1,1],repeat=n)));out['P']=float(np.mean(abs(signs@x/n)>=abs(x.mean())-1e-12));nz=x[x!=0];out['sign_P']=float(stats.binomtest(int((nz>0).sum()),len(nz),.5).pvalue) if len(nz) else 1.;ci=stats.t.interval(.95,n-1,loc=x.mean(),scale=stats.sem(x));out.update(CI_low=float(ci[0]),CI_high=float(ci[1]),LOO_mean_min=float(min(np.delete(x,i).mean() for i in range(n))),LOO_mean_max=float(max(np.delete(x,i).mean() for i in range(n))))
 return out
summ=[]
for key,g in p.groupby(['variant','endpoint','gene']):summ.append(dict(zip(['variant','endpoint','gene'],key),**infer(g.effect),mean_SD=float(g.effect_SD.mean())))
t=pd.DataFrame(summ);t['FDR']=np.nan;t['family']=''
for (variant,role),idx in t.assign(role=np.where(t.gene.eq('CGNL1'),'target','ordinary_comparators')).groupby(['variant','role']).groups.items():
 vals=t.loc[idx,'P'].fillna(1).to_numpy();planned=1 if variant=='matched' and role=='target' else 2 if role=='target' else 6 if variant=='matched' else 12;adj=multipletests(np.r_[vals,np.ones(max(0,planned-len(vals)))],method='fdr_bh')[1][:len(vals)];t.loc[idx,'FDR']=adj;t.loc[idx,'family']=variant+'|'+role+'|endpoints'
 t.loc[t.index.isin(idx)&t.P.isna(),'FDR']=np.nan
t.to_csv(R/'results/all_endpoint_summary.tsv',sep='\t',index=False);t[t.gene.eq('CGNL1')].to_csv(R/'results/CGNL1_endpoint_summary.tsv',sep='\t',index=False)
sp=[]
for (variant,endpoint),g in p.groupby(['variant','endpoint']):
 a=g.pivot(index='patient',columns='gene',values='effect_SD').reindex(columns=targets).dropna()
 if len(a):sp.append(dict(variant=variant,endpoint=endpoint,**infer(a.CGNL1-a[targets[1:]].median(axis=1))))
sp=pd.DataFrame(sp);sp['FDR']=np.nan
for variant,idx in sp.groupby('variant').groups.items():
 vals=sp.loc[idx,'P'].fillna(1).to_numpy();adj=multipletests(np.r_[vals,np.ones(max(0,(1 if variant=='matched' else 2)-len(vals)))],method='fdr_bh')[1][:len(vals)];sp.loc[idx,'FDR']=adj
sp.to_csv(R/'results/target_specificity_summary.tsv',sep='\t',index=False)
assert max(x['max_FWL_error'] for x in checks)<1e-8
done={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'section_effects':len(q),'patient_effects':len(p),'max_FWL_error':max(x['max_FWL_error'] for x in checks),'spec_unchanged':hashlib.sha256((R/'ANALYSIS_SPECIFICATION.json').read_bytes()).hexdigest()==F['sha256'],'no_spot_inferential_P':True}
(R/'ANALYSIS_COMPLETE.json').write_text(json.dumps(done,indent=2),encoding='utf-8');print(t[t.gene.eq('CGNL1')][['variant','endpoint','n','mean','mean_SD','positive','negative','P','FDR','LOO_mean_min','LOO_mean_max']].to_string(index=False));print('SPECIFICITY',sp.to_string(index=False));print(json.dumps(done))
