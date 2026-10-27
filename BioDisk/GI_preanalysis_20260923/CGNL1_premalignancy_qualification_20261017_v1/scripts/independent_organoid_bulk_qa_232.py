
from pathlib import Path
import csv,gzip,json,math,statistics,itertools,datetime
p=Path.cwd();out=p/'results/GSE210993_bulk_230'
with gzip.open(p/'raw/GSE210993_IM.expression_TPM_matrix.103.txt.gz','rt',encoding='utf-8') as h:
 rd=csv.DictReader(h,delimiter=chr(9));rows=[r for r in rd if r['Gene_ID'].split('.')[0]=='ENSG00000128849']
assert len(rows)==1 and rows[0]['Gene_Name']=='CGNL1';raw=rows[0]
meta=json.loads((p/'metadata/GSE210993_ANALYSIS_SAMPLES_227.json').read_text(encoding='utf-8'));lookup={r['sample_name']:r for r in meta}
with (out/'target_samples.tsv').open(encoding='utf-8') as h:saved=list(csv.DictReader(h,delimiter=chr(9)))
assert len(saved)==103 and len({r['sample_name'] for r in saved})==103
error=0.
for r in saved:
 v=float(raw[r['sample_name']]);assert v==float(r['CGNL1_TPM'])
 error=max(error,abs(math.log2(1+v)-float(r['CGNL1_log2TPM1'])))
assert error<1e-12
expected={}
for use_reclassification in [False,True]:
 groups={}
 for r in meta:
  if (r['reclassified'] and not use_reclassification) or r['type'] not in ['Human IM organoid','Human normal gastric organoid']:continue
  groups.setdefault((r['patient_id'],r['type']),[]).append(math.log2(1+float(raw[r['sample_name']])))
 patients=sorted({k[0] for k in groups if (k[0],'Human IM organoid') in groups and (k[0],'Human normal gastric organoid') in groups})
 label='IM_organoid_minus_normal_gastric'+('_author_reclassification_sensitivity' if use_reclassification else '_primary')
 expected[label]={d:statistics.mean(groups[(d,'Human IM organoid')])-statistics.mean(groups[(d,'Human normal gastric organoid')]) for d in patients}
with (out/'matched_specimens.tsv').open(encoding='utf-8') as h:pairrows=list(csv.DictReader(h,delimiter=chr(9)))
grouped={}
for r in pairrows:
 c=r['case_sample'];b=r['reference_sample'];assert c in lookup and b in lookup and lookup[c]['patient_id']==lookup[b]['patient_id']==r['patient_id']
 delta=math.log2(1+float(raw[c]))-math.log2(1+float(raw[b]));assert abs(delta-float(r['difference_log2TPM1']))<1e-12
 grouped.setdefault((r['pair_type'],r['patient_id']),[]).append(delta)
for kind,label in [('long_term','long_term_minus_baseline_IM'),('organoid_minus_tissue','IM_organoid_minus_matched_IM_tissue')]:
 expected[label]={d:statistics.mean(v) for (k,d),v in grouped.items() if k==kind}
def exact_signedrank(d):
 nz=[v for v in d if v!=0];ordered=sorted(enumerate(nz),key=lambda x:abs(x[1]));ranks=[0]*len(nz);position=1
 for absolute,group in itertools.groupby(ordered,key=lambda x:abs(x[1])):
  block=list(group);rank2=2*position+len(block)-1
  for index,value in block:ranks[index]=rank2
  position+=len(block)
 total=sum(ranks);observed=abs(2*sum(r for r,v in zip(ranks,nz) if v>0)-total);distribution={0:1}
 for rank in ranks:
  nxt=dict(distribution)
  for value,count in distribution.items():nxt[value+rank]=nxt.get(value+rank,0)+count
  distribution=nxt
 return sum(count for value,count in distribution.items() if abs(2*value-total)>=observed)/2**len(nz)
with (out/'target_contrasts.tsv').open(encoding='utf-8') as h:contrasts=list(csv.DictReader(h,delimiter=chr(9)))
checks=[]
for r in contrasts:
 delta=list(expected[r['contrast']].values());assert len(delta)==int(r['independent_patients'])
 assert abs(statistics.mean(delta)-float(r['mean_difference_log2TPM1']))<1e-12
 value=exact_signedrank(delta);assert abs(value-float(r['P']))<1e-12
 checks.append({'contrast':r['contrast'],'patients':len(delta),'mean_delta':statistics.mean(delta),'P':value})
secondary=sorted([r for r in contrasts if r['role']=='secondary'],key=lambda r:float(r['P']))
q=[min(2*float(secondary[0]['P']),float(secondary[1]['P'])),float(secondary[1]['P'])]
assert all(abs(float(r['FDR'])-v)<1e-12 for r,v in zip(secondary,q))
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'all103target_values_checked':True,'sample_log_max_error':error,'method':'Independent CSV raw-row lookup; Python scalar log2 and patient means; exact signed-rank distribution via doubled-midrank subset-sum dynamic programming; separate two-test BH arithmetic','checks':checks,'raw_reads_unavailable':'Author processed TPM archive checked; original read alignment not reprocessed'}
(out/'INDEPENDENT_NUMERICAL_QA_232.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
