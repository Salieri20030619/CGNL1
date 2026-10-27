from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,datetime,itertools
import numpy as np,pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path.cwd();S=json.loads((O/'ANALYSIS_SPECIFICATION.json').read_text());covs=S['main_covariates'];targets=S['targets']
note={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'post-primary qualification/robustness, not independent validation','reason':'first matching failed balance in all >=20-pair sections; assess linear-covariate overlap without choosing outcome-favorable matches','method':'unpenalized logistic model for near vs deep using standardized prespecified covariates and macro-component intercepts; overlap weights (1-p) in near and p in deep; outcomes excluded from fitting','qualification':'finite fit, max absolute weighted standardized covariate difference <=0.1, near and deep effective sample size >=20; no relaxed fallback','E1':'weighted near-minus-deep target expression; standardized by original eligible section expression SD','E2':'weighted regression same fully adjusted near x neighbor_mural model as primary; original within-section neighbor standardization','report':'all targets and fixed ordinary comparators; two-endpoint family for CGNL1, separate comparator and specificity families; record failed sections and weight concentration','interpretation':'balance of measured covariate means, not proof of causal exchangeability or cell source identity'}
(O/'OVERLAP_VALIDATION_SPEC.json').write_text(json.dumps(note,indent=2),encoding='utf-8')
d=pd.read_csv(O/'results/spatial_spot_features.tsv.gz',sep='\t');rows=[];qual=[];weightrows=[]
for sample,v in d.groupby('sample'):
 v=v[v.region.eq('Liver macrometastasis tumor')&v.Layer2.eq('Liver_metastatic carcinoma')].copy();v['near']=v.distance.between(50,200).astype(float);v=v[v.near.eq(1)|v.distance.between(500,1500)].dropna(subset=covs+targets+['neighbor_mural_200']).copy();groups=[]
 for k,g in v.groupby('macro_component'):
  if g.near.sum()>=5 and (1-g.near).sum()>=5:groups.append(k)
 v=v[v.macro_component.isin(groups)];a=v.near.to_numpy();base={'sample':sample,'patient':int(v.patient.iloc[0]) if len(v) else int(d[d['sample']==sample].patient.iloc[0]),'near_spots':int(a.sum()),'deep_spots':int(len(a)-a.sum())}
 if min(a.sum(),len(a)-a.sum())<20:qual.append(dict(base,status='insufficient_band_coverage'));continue
 Z=v[covs].to_numpy();sd=Z.std(0);active=sd>1e-10;Z=(Z[:,active]-Z[:,active].mean(0))/sd[active];comp=pd.get_dummies(v.macro_component,drop_first=True,dtype=float).to_numpy();X=np.c_[np.ones(len(v)),Z,comp]
 def fun(b):
  eta=X@b;return float(np.sum(np.logaddexp(0,eta)-a*eta)),X.T@(expit(eta)-a)
 fit=minimize(fun,np.zeros(X.shape[1]),jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':1000});pr=expit(X@fit.x);w=np.where(a==1,1-pr,pr);near=a==1;far=~near
 wn=w[near];wf=w[far];essn=wn.sum()**2/(wn@wn);essf=wf.sum()**2/(wf@wf);balance=np.average(X[near,1:],axis=0,weights=wn)-np.average(X[far,1:],axis=0,weights=wf);maxbal=float(abs(balance).max());grad=float(abs(fit.jac).max());status='qualified' if np.isfinite(w).all() and maxbal<=.1 and min(essn,essf)>=20 and grad<1e-3 else 'insufficient_effective_overlap_or_convergence'
 qual.append(dict(base,status=status,ESS_near=float(essn),ESS_deep=float(essf),max_abs_weighted_SMD=maxbal,max_gradient=grad,optimizer_success=bool(fit.success),optimizer_message=str(fit.message)))
 for barcode,aa,pp,ww in zip(v.barcode,a,pr,w):weightrows.append({'sample':sample,'barcode':barcode,'near':aa,'propensity':pp,'weight':ww,'status':status})
 if status!='qualified':continue
 Y=v[targets].to_numpy();ys=Y.std(0);e1=np.average(Y[near],axis=0,weights=wn)-np.average(Y[far],axis=0,weights=wf)
 nm=v.neighbor_mural_200.to_numpy();nz=(nm-nm.mean())/nm.std();D=np.c_[X,a,nz,a*nz];Dw=D*np.sqrt(w)[:,None];Yw=Y*np.sqrt(w)[:,None];betas=np.linalg.lstsq(Dw,Yw,rcond=None)[0][-1];rank=np.linalg.matrix_rank(Dw);cond=np.linalg.cond(Dw)
 for endpoint,effects in [('E1',e1),('E2',betas)]:
  if endpoint=='E2' and (rank<D.shape[1] or cond>1e8):continue
  for gene,b,s in zip(targets,effects,ys):rows.append(dict(sample=sample,patient=base['patient'],endpoint=endpoint,gene=gene,effect=float(b),effect_SD=float(b/s) if s>1e-10 else np.nan))
pd.DataFrame(qual).to_csv(O/'results/overlap_qualification.tsv',sep='\t',index=False);pd.DataFrame(weightrows).to_csv(O/'results/overlap_weights.tsv.gz',sep='\t',index=False);q=pd.DataFrame(rows);q.to_csv(O/'results/overlap_section_effects.tsv',sep='\t',index=False)
def inf(x):
 x=np.asarray(x);n=len(x);r=dict(n=n,mean=float(x.mean()),positive=int((x>0).sum()),negative=int((x<0).sum()),P=np.nan,CI_low=np.nan,CI_high=np.nan,LOO_min=np.nan,LOO_max=np.nan)
 if n>=2:
  ss=np.array(list(itertools.product([-1,1],repeat=n)));r['P']=float(np.mean(abs(ss@x/n)>=abs(x.mean())-1e-12));ci=stats.t.interval(.95,n-1,loc=x.mean(),scale=stats.sem(x));r.update(CI_low=float(ci[0]),CI_high=float(ci[1]),LOO_min=float(min(np.delete(x,i).mean() for i in range(n))),LOO_max=float(max(np.delete(x,i).mean() for i in range(n))))
 return r
if len(q):
 p=q.groupby(['endpoint','gene','patient'],as_index=False).agg(effect=('effect','mean'),effect_SD=('effect_SD','mean'),sections=('sample','nunique'));p.to_csv(O/'results/overlap_patient_effects.tsv',sep='\t',index=False);summary=[]
 for (endpoint,gene),g in p.groupby(['endpoint','gene']):summary.append(dict(endpoint=endpoint,gene=gene,**inf(g.effect),mean_SD=float(g.effect_SD.mean())))
 t=pd.DataFrame(summary);t['FDR']=np.nan
 for role,idx in t.assign(role=np.where(t.gene=='CGNL1','target','comparator')).groupby('role').groups.items():t.loc[idx,'FDR']=multipletests(t.loc[idx,'P'].fillna(1),method='fdr_bh')[1]
 t.to_csv(O/'results/overlap_summary.tsv',sep='\t',index=False);specific=[]
 for endpoint,g in p.groupby('endpoint'):
  mat=g.pivot(index='patient',columns='gene',values='effect_SD').reindex(columns=targets).dropna();specific.append(dict(endpoint=endpoint,**inf(mat.CGNL1-mat[targets[1:]].median(axis=1))))
 sp=pd.DataFrame(specific);sp['FDR']=multipletests(sp.P.fillna(1),method='fdr_bh')[1];sp.to_csv(O/'results/overlap_specificity.tsv',sep='\t',index=False);print('TARGET',t[t.gene=='CGNL1'].to_string(index=False));print('SPECIFICITY',sp.to_string(index=False))
print('QUALIFICATION',pd.DataFrame(qual).to_string(index=False));(O/'OVERLAP_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'qualified_sections':sum(x['status']=='qualified' for x in qual),'effects':len(q)},indent=2),encoding='utf-8')
