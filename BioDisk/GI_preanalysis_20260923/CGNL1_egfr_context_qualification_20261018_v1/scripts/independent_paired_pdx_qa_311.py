from pathlib import Path
import csv,gzip,json,math,re,datetime,collections
import numpy as np
p=Path.cwd();out=p/'results/paired_pdx_310';a=json.loads((out/'ANALYSIS_COMPLETE_310.json').read_text());spec=json.loads((p/'metadata/STATISTICAL_SPECIFICATION_309.json').read_text());platform={};inside=False;h=None
with (p/'raw/GPL10332_self_297.soft').open(encoding='utf-8') as f:
 for r in csv.reader(f,delimiter='\t'):
  if not r:continue
  if r[0]=='!platform_table_begin':inside=True;continue
  if r[0]=='!platform_table_end':break
  if not inside:continue
  if h is None:h=r;continue
  d=dict(zip(h,r));platform[d['ID']]=d
ids=[];raw=[];header=None;inside=False
with gzip.open(p/'raw/GSE140973_series_matrix.txt.gz','rt',encoding='utf-8') as f:
 for r in csv.reader(f,delimiter='\t'):
  if not r:continue
  if r[0]=='!series_matrix_table_begin':inside=True;continue
  if r[0]=='!series_matrix_table_end':break
  if not inside:continue
  if header is None:header=r[1:];continue
  ids.append(r[0]);raw.append([float(v) for v in r[1:]])
raw=np.asarray(raw);exp=[i for i,k in enumerate(ids) if platform[k]['CONTROL_TYPE']=='FALSE'];sorted_cols=[np.sort(np.log2(raw[exp,j])) for j in range(len(header))];reference=np.vstack(sorted_cols).mean(axis=0);selected={};index={k:i for i,k in enumerate(ids)}
with (out/'normalized_selected_features.tsv').open(encoding='utf-8') as f:stored={r['accession']:r for r in csv.DictReader(f,delimiter='\t')}
maxerr=0.
for probe in spec['probe_manifest']:
 pid=probe['ID'];vals=[]
 for j,gsm in enumerate(header):
  v=math.log2(raw[index[pid],j]);ordered=sorted_cols[j];lo=np.searchsorted(ordered,v,side='left');hi=np.searchsorted(ordered,v,side='right')
  if hi==lo:
   v=float(np.log2(raw[index[pid],j]));lo=np.searchsorted(ordered,v,side='left');hi=np.searchsorted(ordered,v,side='right')
  assert hi>lo;expected=float(reference[lo:hi].mean());err=abs(expected-float(stored[gsm][pid]));maxerr=max(maxerr,err);assert err<1e-10;vals.append(expected)
 selected[pid]=vals
metadata=json.loads((p/'metadata/GSE140973_SAMPLE_METADATA_295.json').read_text());groups=collections.defaultdict(list)
for r in metadata:
 title=r['!Sample_title'][0];m=re.match(r'^BoC([0-9]+)[_ ]+(.*)$',title,re.I);patient='BoC'+m.group(1);token=re.match(r'(5dC|C|K)([0-9]+)',m.group(2),re.I);condition={'5dc':'acute','c':'resistant','k':'untreated'}[token.group(1).lower()];groups[(patient,condition,int(token.group(2)))].append(r['accession'])
assert len(groups)==88 and sum(k[1]=='resistant' for k in groups)==30
symbol_groups=collections.defaultdict(lambda:collections.defaultdict(list))
for r in spec['probe_manifest']:symbol_groups[r['GENE_SYMBOL']][r['NAME']].append(r['ID'])
array={}
for j,gsm in enumerate(header):
 array[gsm]={gene:math.fsum(math.fsum(selected[k][j] for k in v)/len(v) for v in probe_groups.values())/len(probe_groups) for gene,probe_groups in symbol_groups.items()}
patient_values=collections.defaultdict(list)
for (patient,condition,lesion),gsms in groups.items():
 for gene in symbol_groups:patient_values[(patient,condition,gene)].append(math.fsum(array[gsm][gene] for gsm in gsms)/len(gsms))
patient_values={k:math.fsum(v)/len(v) for k,v in patient_values.items()}
with (out/'patient_condition_values.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 for gene in symbol_groups:assert abs(float(r[gene])-patient_values[(r['patient_id'],r['condition'],gene)])<1e-10
patients=sorted({k[0] for k in patient_values});assert len(patients)==10
contrasts={'resistant_minus_acute':('resistant','acute'),'acute_minus_untreated':('acute','untreated'),'resistant_minus_untreated':('resistant','untreated')}
def exact(v):
 v=[round(float(x),12) for x in v];v=[x for x in v if x!=0]
 if not v:return 1.
 ranks={};position=1
 for val in sorted(set(abs(x) for x in v)):
  n=sum(abs(x)==val for x in v);ranks[val]=2*position+n-1;position+=n
 doubled=[ranks[abs(x)] for x in v];total=sum(doubled);observed=abs(2*sum(r for r,x in zip(doubled,v) if x>0)-total);dp={0:1}
 for rank in doubled:
  nxt=dp.copy()
  for score,count in dp.items():nxt[score+rank]=nxt.get(score+rank,0)+count
  dp=nxt
 return sum(count for score,count in dp.items() if abs(2*score-total)>=observed)/(2**len(v))
for res in a['primary_and_secondary']:
 name=res['contrast']
 if res['role']=='lineage_sensitivity':
  v=[patient_values[(pid,'resistant','CGNL1')]-patient_values[(pid,'acute','CGNL1')]-math.fsum(patient_values[(pid,'resistant',g)]-patient_values[(pid,'acute',g)] for g in ['EPCAM','CDH1','KRT8'])/3 for pid in patients]
 else:
  end,start=contrasts[name];v=[patient_values[(pid,end,'CGNL1')]-patient_values[(pid,start,'CGNL1')] for pid in patients]
 assert abs(math.fsum(v)/len(v)-res['mean_difference'])<1e-10 and abs(exact(v)-res['P'])<1e-12
secondary=[r for r in a['primary_and_secondary'] if r['role']=='secondary'];ordered=sorted(secondary,key=lambda r:r['P']);running=1.
for rank in range(len(ordered),0,-1):
 r=ordered[rank-1];running=min(running,r['P']*len(ordered)/rank);assert abs(running-r['FDR'])<1e-12
raw_features={};h=None
with (p/'processed/Agilent_GSM4191604_cp1252_to_utf8_307.txt').open(encoding='utf-8') as f:
 for r in csv.reader(f,delimiter='\t'):
  if not r:continue
  if r[0]=='FEATURES':h=r;continue
  if h and r[0]=='DATA':d=dict(zip(h,r));raw_features[d['FeatureNum']]=d
j=header.index('GSM4191604')
for k,value in zip(ids,raw[:,j]):assert raw_features[k]['ProbeName']==platform[k]['NAME'] and float(raw_features[k]['gProcessedSignal'])==value
neg=[i for i,k in enumerate(ids) if platform[k]['CONTROL_TYPE']=='neg'];qc={}
for j,gsm in enumerate(header):
 vv=sorted(raw[neg,j]);position=.99*(len(vv)-1);lo=int(math.floor(position));expected=vv[lo]+(position-lo)*(vv[lo+1]-vv[lo]);qc[gsm]=expected
with (out/'target_background_QC.tsv').open(encoding='utf-8') as f:
 for r in csv.DictReader(f,delimiter='\t'):assert abs(float(r['negative_control_99pct'])-qc[r['accession']])<1e-8
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'arrays_checked':114,'independent_patients':10,'tumours_checked':len(groups),'resistant_lesions':30,'normalized_feature_values_checked':len(spec['probe_manifest'])*114,'max_normalization_error':maxerr,'raw_pilot_feature_values_checked':len(ids),'DUSP6_acute_decreased_patients':sum(patient_values[(pid,'acute','DUSP6')]<patient_values[(pid,'untreated','DUSP6')] for pid in patients),'methods':'Independent raw gzip parsing, scalar tied-rank quantile reconstruction, title-based lesion/technical aggregation, dynamic-programming exact signed-rank and manual secondary BH; all pilot raw FeatureNum values and negative-control quantiles checked','whole_project_complete':False};(out/'INDEPENDENT_QA_311.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
