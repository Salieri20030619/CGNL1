from pathlib import Path
import json,hashlib,csv,re,datetime,subprocess,warnings,math,sys
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
import statsmodels.api as sm
import statsmodels
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'results/B02_exploratory_scenarios_20261023_v1'
O.mkdir(exist_ok=True)
if (O/'RESULTS.json').exists():raise RuntimeError('Results already exist; do not overwrite a completed exploratory run')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/source_search_593')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
 f=O/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
helper=P/'audit/CGNL1_FOLLOWUP_20261023_v1/followup.py'
ns={};exec(helper.read_text(encoding='utf-8').split('baseline=read')[0],ns);rows=ns['rows']
source_paths=[D/'B02_TableS1.xlsx',D/'B02_protein_NODE.xlsx',D/'B02_phospho_NODE.xlsx',D/'B02_groups_NODE.xlsx',
 P/'audit/CGNL1_FOLLOWUP_20261023_v1/B02_PATIENT_BATCH_PROVENANCE.json',P/'audit/CGNL1_FOLLOWUP_20261023_v1/B02_MODEL_READINESS.json',P/'planning/CGNL1_MASTER_PLAN_20261022_v1/TASKS.json',helper]
before={str(f):sha(f) for f in source_paths};save('BASELINE_HASHES.json',before)
rr=list(rows(D/'B02_TableS1.xlsx','S1A. Clinical info'));h=rr[1][1]
clinical={str(c['A'+str(n)]):{h.get(re.sub(r'\d+$','',k)+'2',k):v for k,v in c.items()} for n,c in rr[2:] if c.get('A'+str(n)) is not None}
provenance=read(P/'audit/CGNL1_FOLLOWUP_20261023_v1/B02_PATIENT_BATCH_PROVENANCE.json')
data=[]
for r in provenance:
 c=clinical[r['case']];bs=sorted(x['batch'] for x in r['channels'])
 assert c['Sex'] in ['Male','Female'] and c['Vascular invasion'] in ['Yes','No']
 data.append(dict(case=r['case'],y=int(c['Vascular invasion']=='Yes'),site=float(r['S283']),protein=float(r['CGNL1_total_protein']),
  age=float(c['Age']),male=int(c['Sex']=='Male'),size_cm=float(c['Tumor_size_diameter (cm)']),
  repeated=int(len(bs)>1),batch_original=bs[0],batch_22=bs[-1]))
df=pd.DataFrame(data);assert len(df)==df.case.nunique()==214 and df.y.sum()==97 and df.repeated.sum()==5
assert np.isfinite(df[['y','site','protein','age','male','size_cm']].to_numpy()).all()

# Verify directly against publisher Table S1; preserve the author-supplied final values.
source_checks=[];source_labels={}
for variable,sheet,rowno,headerrow in [('protein','S1D. proteins expression',3340,3),('site','S1E. 18,347 phosphosites',5520,3)]:
 rr=dict(rows(D/'B02_TableS1.xlsx',sheet,{2,headerrow,rowno}));ids=rr[headerrow];values=rr[rowno];group_cells=rr[2]
 group=None;mapped={};labels={}
 for coord,value in sorted(ids.items(),key=lambda kv:ns['colnum'](kv[0])):
  col=re.sub(r'\d+$','',coord)
  if group_cells.get(col+'2') is not None:group=group_cells[col+'2']
  if isinstance(value,(int,float)):
   case=str(int(value));assert case not in mapped;mapped[case]=values[col+str(rowno)];labels[case]=group
 assert set(mapped)==set(df.case)
 errs=[abs(float(mapped[r.case])-getattr(r,variable)) for r in df.itertuples()]
 assert max(errs)<1e-12
 source_checks.append(dict(variable=variable,sheet=sheet,target_row=rowno,n=214,max_absolute_error=max(errs)))
 source_labels[variable]=labels
repeat_labels=[]
for r in df[df.repeated==1].itertuples():
 expected=int(r.batch_original.split('_')[1]);g=[]
 for variable in ['site','protein']:
  label=source_labels[variable][r.case];assert int(re.search(r'\d+',str(label)).group())==expected;g.append(label)
 repeat_labels.append(dict(case=r.case,original_batch=r.batch_original,publisher_labels=g))
label_discrepancies=[]
for r in df.itertuples():
 for var in ['site','protein']:
  label=source_labels[var][r.case]
  if int(re.search(r'\d+',str(label)).group())!=int(r.batch_original.split('_')[1]):
   label_discrepancies.append(dict(case=r.case,variable=var,publisher_label=label,repository_design=r.batch_original))
save('PAPER_DATA_COMPARISON.json',dict(target_values=source_checks,repeat_batch_labels=repeat_labels,
 other_batch_label_discrepancies=label_discrepancies,handling='Use original NODE channel design for unambiguous batches; preserve publisher header discrepancies explicitly. No target values altered.',
 same_CGNL1_S283_vascular_effect_found_in_maintext_or_figure_legends=False,
 comparator_scope='Exact published sample/value correspondence and repeat-patient batch labels; no matching published CGNL1 effect size available.'))
df['age10']=(df.age-60)/10
df.to_csv(O/'MODEL_INPUT.csv',index=False,encoding='utf-8')
spec=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authorization='User requested running all proposed scenarios and comparison with original paper data/conclusions.',
 scope_change='Explicitly authorized exploratory sensitivity branch despite unresolved source-provenance and mainline gates; prior baseline and formal plan remain unchanged.',
 endpoint='Vascular invasion:Yes=1,No=0',target='S283 log2-relative-abundance coefficient conditional on total CGNL1 and listed covariates',
 family='Four correlated exploratory S283-coefficient tests, one per scenario; BH across four valid tests for supplementary strength assessment, not independent replication.',
 covariates=['total CGNL1 log2 relative abundance','age per10years centered60','male versus female','tumor size linear per cm'],
 input='author observed/NA sheets; 214 unique patients; no new imputation, no dichotomization, no original channel values invented',
 model='unpenalized binomial logistic GLM;MLE;two-sided Wald P and 95% Wald CI',
 scenarios=[dict(id='S0',n=214,batch='none',role='internal no-batch reference'),dict(id='S1',n=209,batch='none',role='exclude five repeated IDs'),
 dict(id='S2',n=214,batch='22-level fixed effect; repeated IDs assigned original batches',role='author-label-aligned exploratory scenario'),
 dict(id='S3',n=214,batch='22-level fixed effect; repeated IDs assigned batch22',role='alternative batch-assignment sensitivity')],
 no_author_effect_to_match=True,selection_criterion='Data fidelity before model results: retain all214 author final values; among batch scenarios use labels displayed for the five repeated IDs in original Table S1. This identifies S2, irrespective of its P or effect direction.',
 selected_for_source_label_fidelity='S2',original_numeric_aggregation_rule_inferred=False,
 result_status='exploratory new association analysis using faithfully reproduced author inputs; not replication of a published CGNL1-specific conclusion',
 model_capacity='Batch scenarios have26 nonintercept parameters for97/117 outcomes; sparse and potentially unstable. Report convergence/rank/influence and do not promote nominal effects to validated evidence.',
 temporal_or_spatial_models='Not fitted: separate timepoint/location observations are unavailable.',
 software=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__),
 fit_results_not_examined_at_freeze=True)
save('SPECIFICATION_BEFORE_FITS.json',spec)
save('SPEC_FREEZE.json',dict(utc=spec['created_utc'],sha256=sha(O/'SPECIFICATION_BEFORE_FITS.json'),input_sha256=sha(O/'MODEL_INPUT.csv')))

all_results=[];coefficients=[];predictions=[];influences=[];designs=[]
for sid in ['S0','S1','S2','S3']:
 sub=df[df.repeated==0].copy() if sid=='S1' else df.copy()
 x=pd.DataFrame({'Intercept':np.ones(len(sub))},index=sub.index)
 for var in ['site','protein','age10','male','size_cm']:x[var]=sub[var]
 if sid in ['S2','S3']:
  col='batch_original' if sid=='S2' else 'batch_22'
  for level in sorted(df[col].unique())[1:]:x['batch_'+level]=(sub[col]==level).astype(float)
 X=x.to_numpy(dtype=float);y=sub.y.to_numpy(dtype=float)
 rank=int(np.linalg.matrix_rank(X));assert rank==X.shape[1]
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always');fit=sm.GLM(y,X,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-11)
 warning_text=[str(w.message) for w in caught]
 beta=float(fit.params[1]);se=float(fit.bse[1]);pvalue=float(2*norm.sf(abs(beta/se)))
 prob=fit.predict();score=np.max(np.abs(X.T@(y-prob)));weights=prob*(1-prob);info=X.T@(weights[:,None]*X)
 valid=bool(fit.converged and np.isfinite(fit.params).all() and np.isfinite(fit.bse).all() and score<1e-5)
 influence=fit.get_influence(observed=True);cooks=influence.cooks_distance[0];hat=influence.hat_matrix_diag
 rec=dict(scenario=sid,n=len(sub),events=int(y.sum()),nonevents=int(len(y)-y.sum()),parameters=X.shape[1],rank=rank,
 minor_class_per_nonintercept_parameter=float(min(y.sum(),len(y)-y.sum())/(X.shape[1]-1)),
 beta=beta,se=se,OR=float(np.exp(beta)),CI95_low=float(np.exp(beta-1.959963984540054*se)),CI95_high=float(np.exp(beta+1.959963984540054*se)),
 P=pvalue if valid else None,fit_valid=valid,converged=bool(fit.converged),iterations=int(fit.fit_history['iteration']),max_score=float(score),
 information_condition_number=float(np.linalg.cond(info)),min_probability=float(prob.min()),max_probability=float(prob.max()),
 max_leverage=float(hat.max()),max_cooks=float(cooks.max()),n_cooks_over_4_over_n=int((cooks>4/len(sub)).sum()),
 warnings=warning_text,interpretation='探索性；批次参数较多' if sid in ['S2','S3'] else '探索性',AIC=float(fit.aic))
 all_results.append(rec)
 for j,name in enumerate(x.columns):coefficients.append(dict(scenario=sid,term=name,beta=float(fit.params[j]),se=float(fit.bse[j]),P=float(fit.pvalues[j])))
 for i,case in enumerate(sub.case):
  predictions.append(dict(scenario=sid,case=case,observed=int(y[i]),probability=float(prob[i])))
  influences.append(dict(scenario=sid,case=case,leverage=float(hat[i]),cooks=float(cooks[i])))
 for row in X:assert np.isfinite(row).all()
 designs.append(dict(scenario=sid,columns=list(x.columns),rows=len(x),rank=rank))
idx=sorted([i for i,r in enumerate(all_results) if r['P'] is not None],key=lambda i:all_results[i]['P']);q=1.
for j in range(len(idx)-1,-1,-1):
 i=idx[j];q=min(q,all_results[i]['P']*len(idx)/(j+1));all_results[i]['FDR_4_scenarios']=float(q)
for r in all_results:r['beta_minus_S0']=r['beta']-all_results[0]['beta']
save('RESULTS.json',all_results);save('DESIGN_DIAGNOSTICS.json',designs)
pd.DataFrame(all_results).to_csv(O/'SCENARIO_RESULTS.csv',index=False,encoding='utf-8')
pd.DataFrame(coefficients).to_csv(O/'ALL_COEFFICIENTS.csv',index=False,encoding='utf-8')
pd.DataFrame(predictions).to_csv(O/'PREDICTIONS.csv',index=False,encoding='utf-8')
pd.DataFrame(influences).to_csv(O/'INFLUENCE.csv',index=False,encoding='utf-8')

# Independently reproduce GLMs in base R; no external packages or auto-selected models.
rscript=r'''args <- commandArgs(trailingOnly=TRUE)
o <- args[1]
d <- read.csv(file.path(o,"MODEL_INPUT.csv"),stringsAsFactors=FALSE)
out <- list()
for(sid in c("S0","S1","S2","S3")) {
 z <- if(sid=="S1") d[d$repeated==0,] else d
 f <- y ~ site + protein + age10 + male + size_cm
 if(sid %in% c("S2","S3")) {
  z$batch <- factor(if(sid=="S2") z$batch_original else z$batch_22,levels=sort(unique(d$batch_original)))
  f <- update(f,. ~ . + batch)
 }
 fit <- glm(f,data=z,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200))
 b <- coef(summary(fit))["site",]
 out[[sid]] <- data.frame(scenario=sid,n=nrow(z),beta=b[1],se=b[2],P=b[4],converged=fit$converged,rank=fit$rank)
}
write.csv(do.call(rbind,out),file.path(o,"R_INDEPENDENT_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(o,"R_SESSION_INFO.txt"),useBytes=TRUE)
'''
(O/'independent_qa.R').write_text(rscript,encoding='utf-8')
proc=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'independent_qa.R'),str(O)],capture_output=True,text=True)
(O/'R_QA_STDOUT.txt').write_text(proc.stdout,encoding='utf-8');(O/'R_QA_STDERR.txt').write_text(proc.stderr,encoding='utf-8');assert proc.returncode==0,proc.stderr
rq=pd.read_csv(O/'R_INDEPENDENT_QA.csv');comparisons=[]
for r in all_results:
 v=rq[rq.scenario==r['scenario']].iloc[0]
 delta={key:abs(float(v[key])-r[key]) for key in ['beta','se','P']}
 assert all(x<1e-6 for x in delta.values()) and bool(v['converged']) and int(v['rank'])==r['rank']
 comparisons.append(dict(scenario=r['scenario'],absolute_differences=delta))
save('INDEPENDENT_NUMERICAL_QA.json',dict(pass_QA=True,method='statsmodels IRLS versus independent base R glm',comparisons=comparisons))
assert all(sha(Path(k))==v for k,v in before.items())
save('ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),models_complete=4,independent_R_QA_pass=True,
 source_target_values_match_paper=428,source_files_unchanged=len(before),selected_source_fidelity_scenario='S2',
 selection_based_on_P_or_direction=False,original_published_CGNL1_effect_reproduced=False,whole_project_complete=False,
 artifacts=[dict(path=str(f),sha256=sha(f)) for f in O.iterdir() if f.is_file() and f.name!='ACCEPTANCE.json']))
print(json.dumps(all_results,ensure_ascii=False));print('QA_PASS',True)
