
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime
import numpy as np,pandas as pd
import statsmodels.api as sm
from scipy import stats
p=Path.cwd()
assert json.loads((p/'metadata/PREPARATION_EXIT_221b.json').read_text())['exit_code']==0
spec=json.loads((p/'STATISTICAL_SPECIFICATION_219.json').read_text());assert (p/'COMMON_GENE_DENOMINATOR_AMENDMENT_221b.json').exists()
models=pd.read_csv(p/'processed/selected_patient_models_221.tsv',sep=chr(9))
assert models.patient_id.is_unique and models.model_id.is_unique
families=spec['families'];drugs=sum(families.values(),[])
aliases={d:d for d in drugs};aliases['5-fluorouracil']='fluorouracil'
choices=[];data_by={};availability=[]
for screen in ['GDSC1','GDSC2']:
 raw=pd.read_csv(p/'raw'/(screen+'_fitted_dose_response_24Jul22.csv'),dtype={'SANGER_MODEL_ID':str,'DRUG_ID':str})
 raw['drug']=raw.DRUG_NAME.str.lower().map(aliases)
 raw=raw[raw.drug.notna()].copy()
 data=raw.merge(models,left_on='SANGER_MODEL_ID',right_on='model_id',how='inner',validate='many_to_one')
 for col in ['MIN_CONC','MAX_CONC','AUC','LN_IC50']:data[col]=pd.to_numeric(data[col],errors='raise')
 for drug in drugs:
  d=data[data.drug==drug].copy()
  family=next(k for k,v in families.items() if drug in v)
  if d.empty:availability.append({'screen':screen,'drug':drug,'family':family,'status':'not_measured_or_no_model_overlap'});continue
  ids=d.groupby('DRUG_ID').patient_id.nunique().reset_index(name='n_patients')
  ids['numeric_id']=ids.DRUG_ID.astype(int);ids=ids.sort_values(['n_patients','numeric_id'],ascending=[False,True])
  drugid=ids.DRUG_ID.iloc[0];sel=d[d.DRUG_ID==drugid]
  ranges=sel.groupby(['MIN_CONC','MAX_CONC']).patient_id.nunique().reset_index(name='n_patients')
  ranges=ranges.sort_values(['n_patients','MIN_CONC','MAX_CONC'],ascending=[False,True,True]);best=ranges.iloc[0]
  sel=sel[(sel.MIN_CONC==best.MIN_CONC)&(sel.MAX_CONC==best.MAX_CONC)]
  resp=sel.groupby('model_id',as_index=False).agg(AUC=('AUC','median'),LN_IC50=('LN_IC50','median'),technical_fitted_curves=('NLME_CURVE_ID','count'))
  fitdata=resp.merge(models,on='model_id',validate='one_to_one');assert fitdata.patient_id.is_unique
  fitdata['screen']=screen;fitdata['drug']=drug;fitdata['family']=family;fitdata['drug_id']=drugid
  data_by[(screen,drug)]=fitdata
  choices.append({'screen':screen,'drug':drug,'family':family,'selected_drug_id':drugid,'selected_MIN_CONC':float(best.MIN_CONC),'selected_MAX_CONC':float(best.MAX_CONC),'matched_patient_models':len(fitdata),'all_ID_coverage':ids[['DRUG_ID','n_patients']].to_dict('records'),'all_selected_ID_ranges':ranges.to_dict('records')})
(p/'metadata/ASSAY_CHOICES_BEFORE_MODEL_FITS_222.json').write_text(json.dumps(choices,indent=2),encoding='utf-8')
results=[];withinlineage=[];overlap=[];data_out=p/'processed/model_fits_222';data_out.mkdir(exist_ok=True)
def fit(d,screen,drug,role,measure='log2CPM1',endpoint='AUC',state_covariates=False,omitted_lineage=None):
 family=next(k for k,v in families.items() if drug in v)
 required=[endpoint,'CGNL1_'+measure]+([g+'_'+measure for g in ['EPCAM','VIM','MKI67']] if state_covariates else [])
 dd=d[np.isfinite(d[required]).all(axis=1)].copy()
 sizes=dd.lineage.value_counts();excluded=sizes[sizes<3].to_dict();dd=dd[dd.lineage.isin(sizes[sizes>=3].index)].copy()
 rec={'screen':screen,'drug':drug,'family':family,'role':role,'expression_measure':measure,'response':endpoint,'n_patient_models':len(dd),'lineages':','.join(sorted(dd.lineage.unique())),'excluded_rare_lineages':json.dumps(excluded),'beta_CGNL1':None,'SE_HC3':None,'P':None,'FDR':None,'CI_low':None,'CI_high':None,'residual_df':None,'status':'not_evaluable'}
 if len(dd)<20:rec['reason']='fewer_than20_patient_models';return rec
 X=pd.DataFrame({'Intercept':np.ones(len(dd)),'CGNL1':dd['CGNL1_'+measure].to_numpy()},index=dd.index)
 X=pd.concat([X,pd.get_dummies(dd[['lineage','RNA_dataset_source']],drop_first=True,dtype=float)],axis=1)
 if state_covariates:
  for g in ['EPCAM','VIM','MKI67']:X[g]=dd[g+'_'+measure].to_numpy()
 y=dd[endpoint].to_numpy()
 if np.linalg.matrix_rank(X.to_numpy())!=X.shape[1] or len(dd)-X.shape[1]<10 or np.ptp(y)==0 or np.ptp(X.CGNL1)==0:
  rec['reason']='rank_variation_or_residual_df_gate';return rec
 fitted=sm.OLS(y,X).fit(cov_type='HC3',use_t=True)
 h=fitted.get_influence().hat_matrix_diag
 if h.max()>=1-1e-10:rec['reason']='unit_leverage';return rec
 ci=fitted.conf_int().loc['CGNL1']
 rec.update(beta_CGNL1=float(fitted.params['CGNL1']),SE_HC3=float(fitted.bse['CGNL1']),P=float(fitted.pvalues['CGNL1']),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),residual_df=float(fitted.df_resid),status='evaluated',reason='',max_leverage=float(h.max()),target_expression_min=float(X.CGNL1.min()),target_expression_max=float(X.CGNL1.max()))
 assert all(np.isfinite(rec[k]) for k in ['beta_CGNL1','SE_HC3','P','CI_low','CI_high'])
 stem=screen+'_'+drug+'_'+role
 X.assign(response=y,model_id=dd.model_id.to_numpy(),patient_id=dd.patient_id.to_numpy()).to_csv(data_out/(stem+'.tsv'),sep=chr(9),index=False,encoding='utf-8')
 return rec
for (screen,drug),d in data_by.items():
 d.to_csv(p/'processed'/(screen+'_'+drug+'_matched_models_222.tsv'),sep=chr(9),index=False,encoding='utf-8')
 for role,measure,endpoint,state in [('primary','log2CPM1','AUC',False),('state_adjusted','log2CPM1','AUC',True),('integer_FPKM_sensitivity','log2FPKM1','AUC',False),('LN_IC50_sensitivity','log2CPM1','LN_IC50',False)]:
  results.append(fit(d,screen,drug,role,measure,endpoint,state))
 for lineage in sorted(d.lineage.unique()):
  results.append(fit(d[d.lineage!=lineage],screen,drug,'leave_out_'+lineage))
  g=d[d.lineage==lineage];rho=stats.spearmanr(g.CGNL1_log2CPM1,g.AUC).statistic if len(g)>=3 and g.CGNL1_log2CPM1.nunique()>1 and g.AUC.nunique()>1 else np.nan
  withinlineage.append({'screen':screen,'drug':drug,'lineage':lineage,'n_models':len(g),'Spearman_rho_descriptive':rho,'P_computed':False})
 for source in sorted(d.RNA_dataset_source.unique()):
  results.append(fit(d[d.RNA_dataset_source==source],screen,drug,'RNA_source_'+('Sanger' if 'Sanger' in source else 'Broad')))
for drug in drugs:
 if ('GDSC1',drug) in data_by and ('GDSC2',drug) in data_by:
  a=data_by[('GDSC1',drug)];b=data_by[('GDSC2',drug)]
  paired=a[['model_id','patient_id','AUC']].merge(b[['model_id','AUC']],on='model_id',suffixes=('_GDSC1','_GDSC2'),validate='one_to_one')
  independent=b[~b.patient_id.isin(a.patient_id)]
  overlap.append({'drug':drug,'GDSC1_models':len(a),'GDSC2_models':len(b),'shared_models':len(paired),'GDSC2_nonoverlap_patient_models':len(independent),'response_cross_assay_rho':float(stats.spearmanr(paired.AUC_GDSC1,paired.AUC_GDSC2).statistic) if len(paired)>=3 else None})
  results.append(fit(independent,'GDSC2',drug,'patient_disjoint_validation'))
r=pd.DataFrame(results)
for (screen,family,role,response,measure),ix in r.groupby(['screen','family','role','response','expression_measure']).groups.items():
 valid=[i for i in ix if pd.notna(r.at[i,'P'])]
 if valid:r.loc[valid,'FDR']=stats.false_discovery_control(r.loc[valid,'P'].to_numpy())
r.to_csv(p/'results/association_models_222.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(withinlineage).to_csv(p/'results/within_lineage_descriptives_222.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(overlap).to_csv(p/'results/cross_assay_overlap_222.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(availability).to_csv(p/'results/unavailable_drugs_222.tsv',sep=chr(9),index=False,encoding='utf-8')
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary':r[r.role=='primary'].replace({np.nan:None}).to_dict('records'),'state_adjusted':r[r.role=='state_adjusted'].replace({np.nan:None}).to_dict('records'),'cross_assay_overlap':overlap,'unavailable':availability,'numeric_QA_pending':True,'causal_mechanism_established':False,'whole_project_complete':False}
(p/'results/ASSOCIATION_ANALYSIS_COMPLETE_222.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(r[r.role.isin(['primary','state_adjusted','patient_disjoint_validation'])][['screen','drug','role','n_patient_models','beta_CGNL1','P','FDR','status']].to_string(index=False),flush=True)
