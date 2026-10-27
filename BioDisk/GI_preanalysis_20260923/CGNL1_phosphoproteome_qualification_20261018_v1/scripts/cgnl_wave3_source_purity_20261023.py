from pathlib import Path
import json,csv,re,hashlib,datetime,collections,subprocess
import numpy as np
import statsmodels.api as sm
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'continuous_exploration_20261023_v1';W=O/'wave3_source_purity';W.mkdir(exist_ok=True)
if (W/'RESULTS.json').exists():raise RuntimeError('Preserve completed wave3')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):
 f=W/name;f.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==x
source_paths=[P/'raw/CGNL1_Q0VF96_380.fasta',P/'metadata/SITE_SEQUENCE_CHECK_380.json',
 P/'raw/PDC_original_384/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Phosphoproteome.phosphosite.tmt11.tsv',
 P/'results/H01_582/SEQUENCE_QA.json',P/'results/B03_recovered_590/SEQUENCE_QA.json',P/'raw/PDAC_phenotype.txt',O/'wave1/patient_contrasts.csv']
before={str(f):sha(f) for f in source_paths}
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Targeted source-localization audit plus limited competing-composition-explanation check, not another endpoint search',
 prior_results_known='Wave1 C01_A and other signals already reviewed;these checks are exploratory follow-up',
 target='C01_A PDAC within-patient delta(S198-S297)',proxies=['WES_purity','WGS_purity'],
 matching='Exact unique patient ID; no guessed aliquot mapping. Reject out-of-range purity units rather than silently rescale.',
 model='Separate OLS contrast~purity*10 for each proxy,HC3 robust covariance,use_t=True; target slope per10 percentagepoints of estimated purity. BH across these two models only.',
 minimum_complete_n=10,spearman='Descriptive rho only;no additional significance testing',
 interpretation='Patient-level computational DNA purity,not biopsy-matched measured cell fractions. Correlation does not prove mediation;null does not exclude normal-tissue composition or cell-type-specific phosphorylation.',
 source_audit='Recheck motif/peptide coordinates and report availability of source localization probabilities. A missing CDAP feature does not prove a harmonized feature is erroneous.',
 no_outcome_endpoint_switch=True,no_imputation=True)
save('SPECIFICATION_BEFORE_RESULTS.json',spec);save('FREEZE.json',dict(sha256=sha(W/'SPECIFICATION_BEFORE_RESULTS.json'),inputs=before))
seq=''.join(l.strip() for l in source_paths[0].read_text().splitlines() if not l.startswith('>'))
motifs=read(source_paths[1])['checks'];motifcheck=[]
for x in motifs:
 if x['source_site'] not in ['S198','S283','S297','S486','S708']:continue
 starts=[m.start() for m in re.finditer('(?='+re.escape(x['motif'])+')',seq)];positions=[s+len(x['motif'])//2+1 for s in starts]
 assert positions==[int(x['source_site'][1:])]
 motifcheck.append(dict(site=x['source_site'],motif=x['motif'],coordinate=positions[0],coordinate_pass=True,localization_probability_in_harmonized_matrix=False))
with source_paths[2].open() as h:
 rd=csv.DictReader(h,delimiter='\t');fields=rd.fieldnames;target=[r for r in rd if r.get('Gene')=='CGNL1']
peptidecheck=[]
for r in target:
 site=r['Phosphosite'].split(':')[-1].upper()
 for pep in r['Peptide'].split(';'):
  starts=[m.start() for m in re.finditer('(?='+re.escape(pep.upper())+')',seq)]
  mod=[[aa.upper()+str(st+j+1) for j,aa in enumerate(pep) if aa.islower()] for st in starts]
  assert mod==[[site]]
  peptidecheck.append(dict(site=site,peptide=pep,unique_sequence_mapping=True,assigned_sites=mod[0],possible_STY_acceptors=sum(aa.upper() in 'STY' for aa in pep)))
export_sites=sorted({r['Phosphosite'].split(':')[-1].upper() for r in target})
save('SITE_MEASUREMENT_AUDIT.json',dict(harmonized_coordinate_checks=motifcheck,PDAC_CDAP_export_sites=export_sites,
 PDAC_CDAP_peptides=peptidecheck,export_annotation_fields=[x for x in fields if not x.endswith('Log Ratio')],
 primary_C01_pair_present_in_alternative_CDAP=all(x in export_sites for x in ['S198','S297']),
 C01='S198/S297 motif coordinates verified;no independent peptide/localization-level confirmation for this pair in downloaded alternative CDAP export. Previously verified S283/S708 sensitivity cannot be generalized to this pair.',
 C02='B03 S486 andS283 both exported as specific modified peptides and prior direct-source reconstruction passed;numeric site localization probability absent from that export.',
 C03='H01 exported S297 assignment maps to a serine-rich peptide;sequence coordinates correct,observation-level localization probabilities not in export. H03 imputation/localization limitations remain.',
 no_data_or_results_corrected=True))
with (O/'wave1/patient_contrasts.csv').open() as h:c={r['case']:float(r['value']) for r in csv.DictReader(h) if r['test']=='C01_A'}
assert len(c)==28
with (P/'raw/PDAC_phenotype.txt').open() as h:ph=list(csv.DictReader(h,delimiter='\t'))
assert len(ph)==len({r['idx'] for r in ph});ph={r['idx']:r for r in ph}
assert set(c)<=set(ph)
results=[];inputs=[];missing=[]
for proxy in spec['proxies']:
 x=[];y=[];ids=[]
 for case,value in c.items():
  raw=ph[case][proxy]
  try:v=float(raw)
  except (ValueError,TypeError):v=np.nan
  if not np.isfinite(v):missing.append(dict(case=case,proxy=proxy,raw_value=raw,reason='source missing'));continue
  assert 0<=v<=1,(proxy,case,v)
  x.append(v*10);y.append(value);ids.append(case);inputs.append(dict(case=case,proxy=proxy,purity=v,value=value))
 n=len(y);assert n>=10 and np.std(x)>0
 X=sm.add_constant(np.asarray(x));fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int()[1]
 results.append(dict(proxy=proxy,n=n,missing=28-n,purity_range=[min(x)/10,max(x)/10],slope_per_10_percentagepoints=float(fit.params[1]),HC3_SE=float(fit.bse[1]),
  CI95=[float(v) for v in ci],P=float(fit.pvalues[1]),descriptive_spearman_rho=float(spearmanr(x,y).statistic),
  r_squared=float(fit.rsquared),max_leverage=float(fit.get_influence().hat_matrix_diag.max()),df_resid=int(fit.df_resid)))
ranked=sorted(results,key=lambda r:r['P']);q=1
for i in range(len(ranked)-1,-1,-1):q=min(q,ranked[i]['P']*len(ranked)/(i+1));ranked[i]['FDR']=q
save('RESULTS.json',results);save('MISSINGNESS.json',missing)
with (W/'inputs.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(inputs[0]));w.writeheader();w.writerows(inputs)
rscript=r'''o<-commandArgs(trailingOnly=TRUE)[1];d<-read.csv(file.path(o,"inputs.csv"),stringsAsFactors=FALSE)
out<-lapply(unique(d$proxy),function(k){z<-d[d$proxy==k,];X<-cbind(1,z$purity*10);y<-z$value;b<-solve(crossprod(X),crossprod(X,y));e<-as.vector(y-X%*%b);bread<-solve(crossprod(X));h<-rowSums((X%*%bread)*X);V<-bread%*%crossprod(X,X*(e/(1-h))^2)%*%bread;se<-sqrt(V[2,2]);data.frame(proxy=k,n=nrow(z),beta=b[2],SE=se,P=2*pt(-abs(b[2]/se),df=nrow(z)-2),rho=cor(z$purity,z$value,method="spearman"))})
write.csv(do.call(rbind,out),file.path(o,"R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
'''
(W/'independent_QA.R').write_text(rscript,encoding='utf-8');proc=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(W/'independent_QA.R'),str(W)],capture_output=True,text=True);assert proc.returncode==0,proc.stderr
with (W/'R_QA.csv').open() as h:
 for r in csv.DictReader(h):
  a=next(x for x in results if x['proxy']==r['proxy']);assert abs(float(r['beta'])-a['slope_per_10_percentagepoints'])<1e-10 and abs(float(r['SE'])-a['HC3_SE'])<1e-10 and abs(float(r['P'])-a['P'])<1e-10
assert all(sha(Path(k))==v for k,v in before.items())
save('ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),R_QA_pass=True,source_files_unchanged=True,tests=2,
 no_independent_validation_claim=True,goal_complete=False,artifacts=[dict(path=str(f),sha256=sha(f)) for f in W.iterdir() if f.is_file() and f.name!='ACCEPTANCE.json']))
print(json.dumps(results));print('ALTERNATIVE_CDAP_PAIR_AVAILABLE',all(x in export_sites for x in ['S198','S297']));print('QA_PASS',True)
