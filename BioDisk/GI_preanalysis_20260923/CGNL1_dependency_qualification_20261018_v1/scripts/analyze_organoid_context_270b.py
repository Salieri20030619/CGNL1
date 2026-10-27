from pathlib import Path
import json,datetime,math
import numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
p=Path.cwd();done=json.loads((p/'metadata/RESOLVED_RNA_COMPLETE_280.json').read_text());assert done['all_requests_pass']
out=p/'results/context_280';out.mkdir(exist_ok=False);base=pd.read_csv(p/'results/organoid_264/target_model_values.tsv',sep='\t',keep_default_na=False);genes={'ENSG00000128849':'CGNL1','ENSG00000143375':'CGN','ENSG00000039068':'CDH1'};records=[];excluded=[]
for r in base.to_dict('records'):
 raw=json.loads((p/'raw/target_rnaseq_resolved_280'/(r['model_id']+'.json')).read_text());vals={};ds=set();sources=set()
 for z in raw['data']:
  a=z['attributes'];assert a['model_id']==r['model_id'];g=genes[a['ensembl_gene_id']];assert g not in vals,('duplicate_gene_record',r['model_id'],g);v=a.get('tpm');vals[g]=float(v) if v is not None else math.nan;ds.add(a['dataset_id']);sources.add(a['data_source'])
 for g in genes.values():r[g+'_log2TPM1']=vals.get(g,math.nan)
 r['RNA_dataset_ids']=','.join(map(str,sorted(ds)));r['RNA_sources']=','.join(sorted(sources));records.append(r)
model=pd.DataFrame(records);model.to_csv(out/'model_source_values.tsv',sep='\t',index=False)
keep=[]
for r in records:
 reason=[]
 if r['tissue'] not in ['Esophagus','Large Intestine','Pancreas','Stomach']:reason.append('outside_prespecified_GI_scope')
 if r['tissue_status']!='Tumour':reason.append('not_annotated_tumour')
 if str(r['qc_pass']).upper()!='TRUE':reason.append('QC_not_passed')
 if not r['patient_id']:reason.append('missing_patient_id')
 if not all(math.isfinite(r[g+'_log2TPM1']) and r[g+'_log2TPM1']>=0 for g in genes.values()):reason.append('incomplete_expression')
 if r['RNA_dataset_ids']!='22' or r['RNA_sources']!='Sanger':reason.append('expression_dataset_requires_separate_qualification')
 if reason:excluded.append({'model_id':r['model_id'],'reason':';'.join(reason)})
 else:keep.append(r)
pd.DataFrame(excluded,columns=['model_id','reason']).to_csv(out/'excluded_models.tsv',sep='\t',index=False)
frame=pd.DataFrame(keep);assert len(frame)>0
patients=[]
for patient,g in frame.groupby('patient_id',sort=True):
 assert g.tissue.nunique()==1,('cross_tissue_patient',patient)
 r={'patient_id':patient,'tissue':g.tissue.iloc[0],'model_ids':';'.join(sorted(g.model_id)),'n_models':len(g),'CGNL1_KO_effect':float(g.CGNL1_fold_change.astype(float).mean())}
 for gene in genes.values():r[gene+'_log2TPM1']=float(g[gene+'_log2TPM1'].mean())
 patients.append(r)
df=pd.DataFrame(patients).set_index('patient_id');df.to_csv(out/'patient_model_values.tsv',sep='\t');assert len(df)>=10
results=[];loo=[]
for gene in ['CGN','CDH1']:
 target=gene+'_log2TPM1';X=pd.concat([pd.DataFrame({'intercept':1.,target:df[target],'CGNL1_log2TPM1':df.CGNL1_log2TPM1},index=df.index),pd.get_dummies(df.tissue,prefix='tissue',drop_first=True,dtype=float)],axis=1);y=df.CGNL1_KO_effect
 assert np.linalg.matrix_rank(X.to_numpy())==X.shape[1];fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc[target];hat=np.sum(X.to_numpy()*(X.to_numpy()@np.linalg.inv(X.T@X)),axis=1)
 for patient in df.index:
  mask=df.index!=patient;sub=sm.OLS(y.loc[mask],X.loc[mask]).fit();loo.append({'target':gene,'omitted_patient':patient,'effect':float(sub.params[target])})
 vals=[r['effect'] for r in loo if r['target']==gene];result={'target':gene,'n_patients':len(df),'n_models':int(df.n_models.sum()),'effect':float(fit.params[target]),'CI_low':float(ci.iloc[0]),'CI_high':float(ci.iloc[1]),'P':float(fit.pvalues[target]),'df_resid':float(fit.df_resid),'max_leverage':float(max(hat)),'condition_number':float(np.linalg.cond(X)),'LOO_effect_min':min(vals),'LOO_effect_max':max(vals),'LOO_same_sign_fraction':sum(np.sign(v)==np.sign(fit.params[target]) for v in vals)/len(vals),'columns':list(X.columns)};results.append(result);pd.concat([y,X],axis=1).to_csv(out/(gene+'_design.tsv'),sep='\t')
fdr=multipletests([r['P'] for r in results],method='fdr_bh')[1]
for r,q in zip(results,fdr):r['FDR']=float(q);r['test_family']='two_prespecified_context_coefficients'
pd.DataFrame(results).to_csv(out/'context_contrasts.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(out/'leave_one_patient_out_effects.tsv',sep='\t',index=False)
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_source_models':len(base),'n_included_models':len(frame),'n_independent_patients':len(df),'by_tissue_patients':df.tissue.value_counts().to_dict(),'results':results,'expression_scale':'API tpm is already log2(TPM+1), used directly','independent_QA_pending':True,'causal_context_validation':False,'whole_project_complete':False};(out/'ANALYSIS_COMPLETE_270.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True)
