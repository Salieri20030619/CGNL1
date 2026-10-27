from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import numpy as np,pandas as pd,json,hashlib
import statsmodels.formula.api as smf
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1');O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
v=pd.read_csv(O/'primary_patient_profiles.tsv',sep='\t');t=pd.read_csv(O/'N1_CGNL1_models.tsv',sep='\t')
assert not v.duplicated(['dataset','definition','patient']).any()
errs=[];fdr=[]
for _,r in t.iterrows():
 u=v[v.dataset.eq(r.dataset)&v.definition.eq(r.definition)&v.cells.ge(r.min_cells)]
 cols=S['models'][r.model]
 fit=smf.ols(r.outcome+' ~ CGNL1 + '+' + '.join(cols),data=u).fit(cov_type='HC3',use_t=True)
 errs.extend([abs(fit.params.CGNL1*u.CGNL1.std()/u[r.outcome].std()-r.beta),abs(fit.pvalues.CGNL1-r.P)])
for _,g in t.groupby('family'):
 p=g.P.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*2/np.arange(1,len(p)+1))[::-1])[::-1]);fdr.extend(abs(q-g.FDR))
assert max(errs)<1e-8 and max(fdr)<1e-10
pr=pd.read_csv(O/'N1_heldout_predictions.tsv',sep='\t');pe=pd.read_csv(O/'N1_predictive_increment.tsv',sep='\t')
keys=['dataset','definition','min_cells','model','outcome']
assert not pr.duplicated(keys+['patient']).any()
cv_errors=[];prediction_errors=[]
for key,g in pr.groupby(keys):
 mask=np.ones(len(pe),bool)
 for k,value in zip(keys,key):mask &= pe[k].eq(value).to_numpy()
 row=pe.loc[mask].iloc[0]
 expected=t
 for k,value in zip(keys,key):expected=expected[expected[k].eq(value)]
 assert len(g)==int(expected.n.iloc[0]),key
 b=(g.baseline-g.observed)**2;a=(g.augmented-g.observed)**2
 cv_errors.append(abs(1-a.mean()/b.mean()-row.improvement))
 if key[1]=='main' and key[2]==20 and key[3]=='M1':
  ds,definition,minimum,model,outcome=key
  u=v[v.dataset.eq(ds)&v.definition.eq(definition)&v.cells.ge(minimum)]
  for _,p in g.iterrows():
   tr=u[u.patient.ne(p.patient)];te=u[u.patient.eq(p.patient)]
   fit=smf.ols(outcome+' ~ CGNL1 + '+' + '.join(S['models'][model]),data=tr).fit()
   prediction_errors.append(abs(float(fit.predict(te).iloc[0])-p.augmented))
assert max(cv_errors)<1e-8 and max(prediction_errors)<1e-8
score=[]
for name,genes in {**S['panels'],**S['outcomes']}.items():
 if all(g in v for g in genes):score.extend(abs(v[genes].mean(1)-v[name]).dropna())
assert max(score)<1e-10
qual=pd.read_csv(O/'primary_cell_coverage.tsv',sep='\t')
for _,r in v.iterrows():
 q=qual[qual.dataset.eq(r.dataset)&qual.definition.eq(r.definition)&qual.patient.eq(r.patient)].iloc[0];assert q.cells==r.cells and q.cells>=20
qa=dict(spec_unchanged=hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256'],
 annotation_amendment_unchanged=hashlib.sha256((R/'ANNOTATION_AMENDMENT.json').read_bytes()).hexdigest()==json.loads((R/'ANNOTATION_AMENDMENT_FREEZE.json').read_text())['sha256'],
 independent_model_error=max(errs),independent_FDR_error=max(fdr),score_reconstruction_error=max(score),prediction_metric_error=max(cv_errors),independent_primary_prediction_error=max(prediction_errors),all_prediction_folds_complete=True,no_duplicate_patients=True,qualified_patients={k:int(g.patient.nunique()) for k,g in v[v.definition.eq('main')].groupby('dataset')},all_pass=True)
assert qa['spec_unchanged'] and qa['annotation_amendment_unchanged']
(R/'NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2));print('QA_COMPLETE',json.dumps(qa),flush=True)
