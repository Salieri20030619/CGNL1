from pathlib import Path
import json,csv,hashlib,datetime,collections,subprocess
import numpy as np
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
ROOT=P/'continuous_exploration_20261023_v1';W=ROOT/'wave2';W.mkdir(exist_ok=True)
if (W/'RESULTS.json').exists():raise RuntimeError('Preserve completed wave2')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
 f=W/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
paths={'P01':('PDAC_phosphosite_357/patient_site_values.tsv','case'),
 'E01':('EC_phosphosite_validation_388/patient_site_values.tsv','patient'),
 'G01':('GASTRIC_validation_415/patient_values.tsv','patient'),
 'G02':('GC2026_validation_418/patient_values.tsv','case'),
 'H01':('H01_582/PATIENT_VALUES.json','case'),
 'H02':('H02_recovered_590/PATIENT_VALUES.json','case'),
 'H03':('HCC_exploratory_572/patient_site_values572.json','patient')}
before={str(P/'results'/rel):sha(P/'results'/rel) for rel,k in paths.values()}
tests=[]
for c in ['G01','G02']:
 for a,b in [('S198','S297'),('S283','S708')]:tests.append(dict(test='C01_'+c+'_'+a+'_'+b,family='C01_cross',cohort=c,a=a,b=b))
for c in ['P01','E01','H01','H02','H03']:tests.append(dict(test='C02_'+c+'_S486_S283',family='C02_cross',cohort=c,a='S486',b='S283'))
for c in ['P01','G01','G02']:tests.append(dict(test='C03_'+c+'_S283_S297',family='C03_cross',cohort=c,a='S283',b='S297'))
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),purpose='Evaluate the same candidate site contrasts across all existing eligible paired data; do not choose a comparison based on favorable P.',
 prior_results='Wave1 already examined. Wave2 is explicitly exploratory follow-up and not blinded independent confirmation.',
 tests=tests,method='Same two-sided exact signed rank and bootstrap mean CI as wave1; preserve original n<10 descriptive convention; no site substitution or merged-site splitting.',
 FDR='BH separately among valid new tests in C01_cross,C02_cross,C03_cross; preserve wave1 FDR separately; correlated reuse of datasets is not independent replication.',
 no_between_cohort_inferential_test=True,no_cancer_specificity_from_P_difference=True,
 independent_PDAC_requalification='FudanTableS3D only aggregate target siteFC/P; Britton has S283 and multiply phosphorylated S283/S293,S296/S297,S486/S491. Neither currently supplies the two exact paired primary endpoints. No substitutions or newP from summaries.',
 normalization='Use source-locked patient-level measurements; site-minus-site cancels a shared denominator only where present. Do not pool effects across unverified measurement scales.')
save('SPEC_BEFORE_RESULTS.json',spec);save('FREEZE.json',dict(sha256=sha(W/'SPEC_BEFORE_RESULTS.json'),inputs=before))
def signed(a):
 y=[float(x) for x in a if x!=0]
 if not y:return None
 rr=[int(round(x*2)) for x in rankdata(np.abs(y))];obs=sum(r for r,x in zip(rr,y) if x>0);cut=min(obs,sum(rr)-obs)
 dp=[0]*(cut+1);dp[0]=1
 for r in rr:
  for j in range(cut,r-1,-1):dp[j]+=dp[j-r]
 return min(1.,2*sum(dp)/(2**len(y)))
values={}
for c,(rel,key) in paths.items():
 f=P/'results'/rel
 if f.suffix=='.tsv':
  with f.open(encoding='utf-8') as h:rr=list(csv.DictReader(h,delimiter='\t'))
 else:rr=read(f)
 m=collections.defaultdict(dict)
 for r in rr:
  site=r.get('site',r.get('feature'));site=site.split('|')[2] if '|' in site else site
  # Combined-site features remain combined and cannot satisfy an exact S297 endpoint.
  v=r.get('adjusted_delta',r.get('adjusted'));case=r[key];assert case not in m[site];m[site][case]=float(v)
 values[c]=m
results=[];patient_rows=[];rng=np.random.default_rng(20261024)
for t in tests:
 m=values[t['cohort']];ids=sorted(set(m[t['a']])&set(m[t['b']]));x=np.array([m[t['a']][c]-m[t['b']][c] for c in ids]);n=len(x)
 for c,v in zip(ids,x):patient_rows.append(dict(test=t['test'],cohort=t['cohort'],case=c,value=float(v)))
 rec=dict(t,n=n,mean=float(x.mean()) if n else None,median=float(np.median(x)) if n else None,
  up=int(sum(x>0)),down=int(sum(x<0)),P=signed(x) if n>=10 else None,
  status='testable_exploratory' if n>=10 else ('descriptive_low_coverage' if n else 'no_common_exact_sites'))
 if n>=10:
  boot=x[rng.integers(0,n,size=(10000,n))].mean(axis=1);rec['bootstrap_mean_CI95']=[float(v) for v in np.quantile(boot,[.025,.975])]
 results.append(rec)
for fam in ['C01_cross','C02_cross','C03_cross']:
 a=sorted([r for r in results if r['family']==fam and r['P'] is not None],key=lambda r:r['P']);q=1.
 for i in range(len(a)-1,-1,-1):q=min(q,a[i]['P']*len(a)/(i+1));a[i]['FDR']=q
save('RESULTS.json',results)
with (W/'patient_contrasts.csv').open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(patient_rows[0]));w.writeheader();w.writerows(patient_rows)
rscript=r'''o<-commandArgs(trailingOnly=TRUE)[1];d<-read.csv(file.path(o,"patient_contrasts.csv"),stringsAsFactors=FALSE)
out<-lapply(unique(d$test),function(k){x<-d$value[d$test==k];n<-length(x);p<-if(n>=10) wilcox.test(x,exact=TRUE,correct=FALSE)$p.value else NA_real_;data.frame(test=k,n=n,mean=mean(x),median=median(x),P=p)})
write.csv(do.call(rbind,out),file.path(o,"R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
'''
(W/'independent_QA.R').write_text(rscript,encoding='utf-8');rr=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(W/'independent_QA.R'),str(W)],capture_output=True,text=True)
(W/'R_stderr.txt').write_text(rr.stderr,encoding='utf-8');assert rr.returncode==0,rr.stderr
with (W/'R_QA.csv').open() as h:qa=list(csv.DictReader(h))
for q in qa:
 r=next(x for x in results if x['test']==q['test']);assert abs(float(q['mean'])-r['mean'])<1e-10
 if r['P'] is not None:assert abs(float(q['P'])-r['P'])<1e-10
assert all(sha(Path(k))==v for k,v in before.items())
save('ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),R_QA_pass=True,
 valid_new_tests=sum(r['P'] is not None for r in results),sources_unchanged=True,goal_complete=False,
 artifacts=[dict(path=str(f),sha256=sha(f)) for f in W.iterdir() if f.is_file() and f.name!='ACCEPTANCE.json']))
save('INDEPENDENT_PDAC_SOURCE_REQUALIFICATION.json',dict(
 Fudan=dict(file=str(P/'raw/Fudan_quantification_TableS3_363.xlsx'),sha256=sha(P/'raw/Fudan_quantification_TableS3_363.xlsx'),status='not_evaluable_from_published_supplement',reason='Only aggregate S283 FC/P; no patient-level S198/S297 or S283/S708 joint site measurements in TableS3.'),
 Britton=dict(source=str(P/'metadata/BRITTON_TARGET_EXTRACTION_396.json'),status='not_evaluable_for_locked_exact_primary_pairs',observed_features=['S283','S283;S293','S296;S297','S486;S491'],reason='S198 and S708 absent from extracted target features; multiply phosphorylated peptides not substituted as isolated S297/S486.'),
 changed_endpoint='Relative site contrasts do not need protein denominator, but still need matching individual site measurements and exact localization; absent inputs were not fabricated.'))
print(json.dumps(results,ensure_ascii=False));print('QA_PASS',True)
