from pathlib import Path
import json,csv,hashlib,datetime,collections,statistics,re,subprocess,sys
import numpy as np
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'continuous_exploration_20261023_v1';O.mkdir(exist_ok=True)
W=O/'wave1';W.mkdir(exist_ok=True)
if (W/'RESULTS.json').exists():raise RuntimeError('Preserve completed results; use another version for changes')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,d):
 f=W/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
files={
 'P01':P/'results/PDAC_phosphosite_357/patient_site_values.tsv',
 'H01':P/'results/H01_582/PATIENT_VALUES.json',
 'H02':P/'results/H02_recovered_590/PATIENT_VALUES.json',
 'H03':P/'results/HCC_exploratory_572/patient_site_values572.json',
 'B03':P/'results/B03_recovered_590/PATIENT_VALUES.json',
 'PDAC_alt':P/'results/PDC_processing_sensitivity_399/patient_values.tsv'}
before={str(f):sha(f) for f in files.values()}
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),goal='C01-C04 sustained exploration authorized by user',
 prior_knowledge='Candidate sites selected after review of existing marginal site effects. These are exploratory new contrasts, not blinded confirmatory hypotheses.',
 primary_tests=[dict(id='C01_A',cohort='P01',contrast='S198 minus S297',expected_n=28,family='C01'),dict(id='C01_B',cohort='P01',contrast='S283 minus S708',expected_n=21,family='C01'),
 dict(id='C02_A',cohort='B03',contrast='S486 protein-adjusted tumor-normal; recomputation/robustness of known effect',expected_n=60,family='C02'),
 dict(id='C02_B',cohort='B03',contrast='S486 minus S283',expected_n=60,family='C02'),
 dict(id='C03_H01',cohort='H01',contrast='S283 minus S297',expected_n=34,family='C03'),dict(id='C03_H02',cohort='H02',contrast='S283 minus S297',expected_n=17,family='C03'),dict(id='C03_H03',cohort='H03',contrast='S283 minus S297',expected_n=27,family='C03')],
 endpoint='Within each patient: (siteA_T-siteA_N)-(siteB_T-siteB_N). Use only jointly observed eligible patient-level values. For C02_A retain original adjusted endpoint.',
 method='Two-sided exact signed-rank using doubled average ranks; zero values excluded; all-zero P unavailable. BH separately within C01(2), C02(2), C03(3). RawP<.05 denotes nominal signal. No pooling of HCC patients.',
 uncertainty='Patient bootstrap percentile95%CI of mean,10000 resamples,seed20261023. Descriptive inference limitations retained.',
 sensitivity=['compare centered and uncentered contrast at identical patients','leave-one-patient-out mean and exactP ranges, not used for exclusion','PDAC alternative CDAP processing for S283-S708 on common patients,not independent validation'],
 C02_source_audit='Reconstruct S486 site and total protein components from source target rows, saved verified offsets, and original case/sample map; assess coverage60/104 and sequence assignment.',
 C03_boundary='Compare within-cohort results descriptively; do not test cross-cohort disease-stage interactions while overlap and source comparability remain unresolved.',
 C04_gate='Only source qualification until gene/site identity, sample-level site+protein values and biological replication verified. Species mapping required; RNA/totalproteome-only sources excluded from phosphorylation inference.',
 new_missing_value_imputation=False,posthoc_outlier_removal=False,new_results_unseen=True)
if (W/'SPECIFICATION_BEFORE_RESULTS.json').exists():
 oldspec=read(W/'SPECIFICATION_BEFORE_RESULTS.json')
 assert {k:v for k,v in oldspec.items() if k!='utc'}=={k:v for k,v in spec.items() if k!='utc'}
 spec=oldspec
 assert read(W/'SPEC_FREEZE.json')['input_hashes']==before
else:
 save('SPECIFICATION_BEFORE_RESULTS.json',spec);save('SPEC_FREEZE.json',dict(utc=spec['utc'],sha256=sha(W/'SPECIFICATION_BEFORE_RESULTS.json'),input_hashes=before))
def signed(vals):
 y=[float(x) for x in vals if x!=0]
 if not y:return None
 ranks=[int(round(2*x)) for x in rankdata(np.abs(y))];obs=sum(r for r,v in zip(ranks,y) if v>0);cut=min(obs,sum(ranks)-obs)
 ways=[0]*(cut+1);ways[0]=1
 for r in ranks:
  for j in range(cut,r-1,-1):ways[j]+=ways[j-r]
 return min(1.,2*sum(ways)/(2**len(y)))
patients={}
for cohort,f in files.items():
 if cohort=='PDAC_alt':continue
 if f.suffix=='.tsv':
  with f.open(encoding='utf-8') as h:rr=list(csv.DictReader(h,delimiter='\t'))
 else:rr=read(f)
 mapped=collections.defaultdict(dict)
 for r in rr:
  site=r['site'].split('|')[2] if '|' in r['site'] else r['site'];case=r.get('case',r.get('patient'))
  assert case not in mapped[site]
  if cohort=='P01':v=float(r['adjusted_delta']);u=float(r['site_delta'])-float(r['protein_delta'])
  elif cohort=='H01':v=r['adjusted'];u=r['unscaled_adjusted']
  elif cohort=='H03':v=r['adjusted_delta'];u=r['uncentered_adjusted_delta']
  else:v=r['adjusted'];u=r['uncentered']
  mapped[site][case]=dict(adjusted=v,uncentered=u,original=r)
 patients[cohort]=mapped
input_rows=[];results=[];rng=np.random.default_rng(20261023)
for test in spec['primary_tests']:
 cohort=test['cohort'];mm=patients[cohort]
 pair=test['contrast'].split(' minus ') if ' minus ' in test['contrast'] else ['S486']
 ids=sorted(set(mm[pair[0]]) & set(mm[pair[1]])) if len(pair)==2 else sorted(mm[pair[0]])
 assert len(ids)==test['expected_n']
 vals=[];unc=[]
 for case in ids:
  va=mm[pair[0]][case];vb=mm[pair[1]][case] if len(pair)==2 else dict(adjusted=0,uncentered=0)
  v=va['adjusted']-vb['adjusted'];u=va['uncentered']-vb['uncentered']
  if cohort=='P01':
   ra=va['original'];rb=vb['original']
   assert abs(float(ra['protein_T'])-float(rb['protein_T']))<1e-12 and abs(float(ra['protein_N'])-float(rb['protein_N']))<1e-12
   u=(float(ra['site_T'])-float(ra['site_N']))-(float(rb['site_T'])-float(rb['site_N']))
  vals.append(v);unc.append(u);input_rows.append(dict(test=test['id'],family=test['family'],cohort=cohort,case=case,value=v,alternative=u))
 a=np.asarray(vals);b=np.asarray(unc);boot=np.mean(a[rng.integers(0,len(a),size=(10000,len(a)))],axis=1)
 loo=[np.delete(a,i) for i in range(len(a))];lp=[signed(v) for v in loo];means=[float(np.mean(v)) for v in loo]
 maxdiff=float(np.max(np.abs(a-b)))
 results.append(dict(test=test['id'],family=test['family'],cohort=cohort,contrast=test['contrast'],n=len(a),
  mean=float(np.mean(a)),median=float(np.median(a)),up=int((a>0).sum()),down=int((a<0).sum()),zeros=int((a==0).sum()),
  bootstrap_mean_CI95=[float(x) for x in np.quantile(boot,[.025,.975])],P=signed(a),
  alternative_mean=float(np.mean(b)),alternative_P=signed(b),max_centering_difference=maxdiff,
  normalization_invariance_expected=len(pair)==2,normalization_invariance_pass=maxdiff<1e-9 if len(pair)==2 else None,
  LOO_mean_range=[min(means),max(means)],LOO_P_range=[min(x for x in lp if x is not None),max(x for x in lp if x is not None)],
  LOO_mean_direction_constant=all(np.sign(x)==np.sign(np.mean(a)) for x in means),
  ties_in_absolute_nonzero_ranks=bool(len(set(abs(x) for x in a if x!=0))<sum(a!=0))))
for family in ['C01','C02','C03']:
 fam=sorted([r for r in results if r['family']==family and r['P'] is not None],key=lambda r:r['P']);q=1.
 for i in range(len(fam)-1,-1,-1):q=min(q,fam[i]['P']*len(fam)/(i+1));fam[i]['FDR']=q
save('RESULTS.json',results)
with (W/'patient_contrasts.csv').open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(input_rows[0]));w.writeheader();w.writerows(input_rows)

# Same-patient CDAP reprocessing check: limited to existing S283/S708, not a new cohort.
with files['PDAC_alt'].open(encoding='utf-8') as h:altrows=list(csv.DictReader(h,delimiter='\t'))
alt=[]
for method in ['unshared','all_peptides']:
 m=collections.defaultdict(dict)
 for r in altrows:
  if r['method']==method:m[r['site']][r['case']]=float(r['adjusted_delta'])
 ids=sorted(set(m['S283'])&set(m['S708'])&set(patients['P01']['S283'])&set(patients['P01']['S708']))
 x=[m['S283'][c]-m['S708'][c] for c in ids]
 base=[patients['P01']['S283'][c]['adjusted']-patients['P01']['S708'][c]['adjusted'] for c in ids]
 alt.append(dict(method=method,n=len(ids),mean=float(np.mean(x)),P=signed(x),base_mean_same_patients=float(np.mean(base)),base_P_same_patients=signed(base),
  individual_sign_agreement=int(sum(np.sign(a)==np.sign(b) for a,b in zip(x,base))),not_independent_validation=True))
save('PDAC_PROCESSING_SENSITIVITY.json',alt)

# B03: direct source target reconstruction, retaining source units and both normalization choices.
source_manifest=read(P/'results/B03_recovered_590/SOURCE_MANIFEST.json')
offsets=read(P/'results/B03_recovered_590/NORMALIZATION.json');pairing=read(P/'results/B03_recovered_590/PAIRING.json')['eligible']
target={};heads={}
for kind,source in zip(['protein','phosphosite'],source_manifest):
 f=Path(source['path']);assert sha(f)==source['sha256'];before[str(f)]=source['sha256']
 with f.open(encoding='utf-8') as h:
  rd=csv.DictReader(h,delimiter='\t');heads[kind]=[x for x in rd.fieldnames if x.startswith('Log ')];target[kind]=[x for x in rd if x.get('Gene')=='CGNL1']
assert len(target['protein'])==1
pr=target['protein'][0];sr=next(x for x in target['phosphosite'] if x['Phosphosite'].split(':')[-1].lower()=='s486')
seq=''.join(x.strip() for x in (P/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not x.startswith('>'))
mapping=[]
for peptide in sr['Peptide'].split(';'):
 starts=[m.start() for m in re.finditer('(?='+re.escape(peptide.upper())+')',seq)]
 sitecoords=[[aa.upper()+str(st+j+1) for j,aa in enumerate(peptide) if aa.islower()] for st in starts]
 mapping.append(dict(peptide=peptide,coordinates=sitecoords,unique_single_S486=sitecoords==[['S486']]))
assert all(x['unique_single_S486'] for x in mapping)
def sample(row,kind,label,center):
 vals=[]
 for c in heads[kind]:
  if c[4:].split(':')[0]!=label:continue
  try:v=float(row[c])
  except (ValueError,TypeError):continue
  if np.isfinite(v):vals.append(v-(offsets[kind][c] if center else 0))
 return float(np.median(vals)) if vals else None
components=[]
for pa in pairing:
 if pa['case'] not in patients['B03']['S486']:continue
 out=dict(case=pa['case'])
 for center,name in [(True,'centered'),(False,'source')]:
  st,sn=[sample(sr,'phosphosite',pa[k],center) for k in ['tumor','normal']]
  pt,pn=[sample(pr,'protein',pa[k],center) for k in ['tumor','normal']]
  value=st-sn-pt+pn;expected=patients['B03']['S486'][pa['case']]['adjusted' if center else 'uncentered'];assert abs(value-expected)<1e-10
  out[name]=dict(site_delta=st-sn,protein_delta=pt-pn,adjusted=value)
 components.append(out)
assert len(components)==60
save('C02_SOURCE_COMPONENTS.json',dict(n=60,eligible_patients=len(pairing),sequence_checks=mapping,components=components,
 source_limit='Single-modification peptide sequence mapping verified; export lacks observation-level localization probabilities. This is not occupancy.',
 coverage_inference='Complete-case result, not evidence for all104 eligible patients. No new imputation.'))

rscript=r'''args<-commandArgs(trailingOnly=TRUE);o<-args[1];d<-read.csv(file.path(o,"patient_contrasts.csv"),stringsAsFactors=FALSE)
exactp<-function(x) {
 x<-x[x!=0];if(!length(x))return(NA_real_)
 if(!anyDuplicated(abs(x)))return(wilcox.test(x,mu=0,alternative="two.sided",exact=TRUE,correct=FALSE)$p.value)
 ranks<-as.integer(round(2*rank(abs(x))));obs<-sum(ranks[x>0]);cut<-min(obs,sum(ranks)-obs)
 p<-c(1,rep(0,cut))
 for(r in ranks){old<-p;p<-old/2;if(r<=cut)p[(r+1):(cut+1)]<-p[(r+1):(cut+1)]+old[1:(cut-r+1)]/2}
 min(1,2*sum(p))
}
out<-lapply(unique(d$test),function(k){z<-d[d$test==k,];data.frame(test=k,n=nrow(z),mean=mean(z$value),median=median(z$value),P=exactp(z$value),alternative_P=exactp(z$alternative))})
write.csv(do.call(rbind,out),file.path(o,"independent_R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(o,"R_session.txt"),useBytes=TRUE)
'''
(W/'independent_QA.R').write_text(rscript,encoding='utf-8')
run=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(W/'independent_QA.R'),str(W)],capture_output=True,text=True)
(W/'R_stderr.txt').write_text(run.stderr,encoding='utf-8');assert run.returncode==0,run.stderr
with (W/'independent_R_QA.csv').open() as h:rr=list(csv.DictReader(h))
qa=[]
for row in rr:
 r=next(x for x in results if x['test']==row['test']);diff={k:abs(float(row[k])-r[k]) for k in ['mean','median','P','alternative_P']};assert max(diff.values())<1e-10
 qa.append(dict(test=r['test'],absolute_difference=diff))
assert all(sha(Path(k))==v for k,v in before.items())
save('INDEPENDENT_QA.json',dict(all_pass=True,checks=qa,C02_source_reconstruction_pass=True,input_hashes_unchanged=before))
save('ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),tests=7,QA_pass=True,
 C04='source qualification ongoing',goal_complete=False,whole_project_complete=False,
 artifacts=[dict(path=str(f),sha256=sha(f)) for f in W.iterdir() if f.is_file() and f.name!='ACCEPTANCE.json']))
print(json.dumps(results,ensure_ascii=False));print('PROCESSING_SENSITIVITY',json.dumps(alt));print('QA_PASS',True)
