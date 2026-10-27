from pathlib import Path
import json,hashlib,datetime,csv,statistics,random,math
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024';W=O/'wave4_phosphoform';W.mkdir(exist_ok=True)
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
inputs=[D/'rollup_audit/TARGET_SOURCE_ROWS.json',D/'lineage/SAMPLE_MAPPING.json']
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),id='C01_PHOSPHOFORM_01',prior_knowledge='Exploratory refinement after known marginal site effects and sourceformaudit,notblindedconfirmation',endpoint='Withinpatient tumor-normal difference in log2 doubleS198S202 minus monoS202,usingUMichmultisiteexports',inclusion='Allcaseswithuniqueofficial mapping and all4values;no newimputation or data-dependentexclusion',test='two-sided exactsignedrank;zerosremoved;doubledaverageranksfor ties',family='one newlylockedexploratory contrast;originalfamilies unchanged;no newFDRsignificanceclaim',sensitivity='BCM S198minusS202 at same patients;LOOmeans;no selectionbyP',uncertainty='10000patientbootstrap mean percentile95CI seed20261024',sources={str(f):sha(f) for f in inputs},hypothesis_selected_before_new_contrast_results=True)
sf=W/'SPEC_BEFORE_RESULTS.json'
if not sf.exists():save(sf,spec)
else:assert json.loads(sf.read_text())['sources']==spec['sources']
if __import__('sys').argv[-1]=='freeze':print('Frozen',sha(sf));raise SystemExit
rows=json.loads(inputs[0].read_text())['multi-site'];mapping=json.loads(inputs[1].read_text());maps={(x['case'],x['tissue']):x['source_column'] for x in mapping}
double=next(x for x in rows if x['Index'].endswith('_2_2_S198S202'));mono=next(x for x in rows if x['Index'].endswith('_1_1_S202'))
bcm={}
for tissue in ['Tumor','Normal']:
 with (D/'BCM'/f'PDAC_phospho_site_abundance_log2_reference_intensity_normalized_{tissue}.txt').open() as h:
  for row in csv.DictReader(h,delimiter='\t'):
   if row['idx'].startswith('ENSG00000128849.') and row['idx'].split('|')[2] in ['S198','S202']:bcm[(tissue,row['idx'].split('|')[2])]=row
records=[];excluded=[]
cases=sorted(set(c for c,t in maps if t=='Tumor')&set(c for c,t in maps if t=='Normal'))
for case in cases:
 tc=maps[(case,'Tumor')];nc=maps[(case,'Normal')];vs=[double[tc],double[nc],mono[tc],mono[nc]]
 if 'NA' in vs:excluded.append(dict(case=case,components_present=[x!='NA' for x in vs]));continue
 a,b,c,d=map(float,vs);alt=float(bcm['Tumor','S198'][case])-float(bcm['Normal','S198'][case])-float(bcm['Tumor','S202'][case])+float(bcm['Normal','S202'][case])
 records.append(dict(test='C01_PHOSPHOFORM_01',case=case,double_T=a,double_N=b,mono_T=c,mono_N=d,value=a-b-c+d,alternative=alt))
with (W/'patient_contrasts.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
def exact(x):
 x=[v for v in x if v!=0];order=sorted(range(len(x)),key=lambda i:abs(x[i]));ranks=[0]*len(x);i=0
 while i<len(x):
  j=i+1
  while j<len(x) and abs(x[order[j]])==abs(x[order[i]]):j+=1
  for k in range(i,j):ranks[order[k]]=i+j+1
  i=j
 obs=sum(r for r,v in zip(ranks,x) if v>0);cut=min(obs,sum(ranks)-obs);counts=[1]+[0]*cut
 for r in ranks:
  for k in range(cut,r-1,-1):counts[k]+=counts[k-r]
 return min(1,2*sum(counts)/(2**len(x))) if x else None
x=[r['value'] for r in records];n=len(x);rng=random.Random(20261024);boot=sorted(statistics.mean(rng.choices(x,k=n)) for _ in range(10000))
def q(p):
 t=(len(boot)-1)*p;i=int(t);return boot[i]+(boot[min(i+1,len(boot)-1)]-boot[i])*(t-i)
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),test='C01_PHOSPHOFORM_01',paired_available=len(cases),n=n,excluded=len(excluded),mean=statistics.mean(x),median=statistics.median(x),down=sum(v<0 for v in x),up=sum(v>0 for v in x),P=exact(x) if n>=10 else None,FDR=None,FDR_reason='single newlylockedexploratorycontrast;retainoriginalfamilies',bootstrap_mean_CI95=[q(.025),q(.975)],LOO_mean_range=[min((sum(x)-v)/(n-1) for v in x),max((sum(x)-v)/(n-1) for v in x)],alternative_BCM_mean=statistics.mean(r['alternative'] for r in records),alternative_P=exact([r['alternative'] for r in records]) if n>=10 else None,not_occupancy=True,not_independent_validation=True)
save(W/'RESULTS.json',result);save(W/'EXCLUSIONS.json',excluded)
(W/'independent_QA.R').write_text((O/'wave1/independent_QA.R').read_text())
print(json.dumps(result))
