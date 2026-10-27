from pathlib import Path
import gzip,csv,json,math,statistics,datetime,collections
O=Path.cwd()/'mechanism_investigation_v2';mapping=json.loads((O/'metadata/GRCh37_CGNL1_expanded.json').read_text());ids={r['id'] for r in mapping['Transcript']};pc={r['id'] for r in mapping['Transcript'] if r['biotype']=='protein_coding'}
with (O/'results/GSE237876_sample_values.tsv').open(encoding='utf-8') as f:samples=list(csv.DictReader(f,delimiter='\t'))
with (O/'results/GSE237876_paired_deltas.tsv').open(encoding='utf-8') as f:deltas=list(csv.DictReader(f,delimiter='\t'))
plan={r['accession']:r for r in json.loads((O/'metadata/ROUTE_INPUT_PLAN_43.json').read_text())};groups=collections.defaultdict(list);error=0.
for row in samples:
 p=O/'raw/GSE237876'/plan[row['sample']]['files'][0].rsplit('/',1)[1];v={}
 with gzip.open(p,'rt',encoding='utf-8') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['Name'] in ids:assert r['Name'] not in v;v[r['Name']]=(float(r['TPM']),float(r['NumReads']))
 assert set(v)==ids
 total=sum(t[0] for t in v.values());coding=sum(v[k][0] for k in pc);log=math.log2(total+1);clog=math.log2(coding+1);error=max(error,abs(total-float(row['TPM'])),abs(coding-float(row['coding_TPM'])),abs(log-float(row['log2TPM1'])));groups[(row['patient'],row['site'])].append((log,clog))
means={k:(statistics.mean(x[0] for x in v),statistics.mean(x[1] for x in v)) for k,v in groups.items()}
for r in deltas:
 g=means[(r['patient'],r['site'])];p=means[(r['patient'],'primary')];error=max(error,abs(g[0]-p[0]-float(r['delta_log2TPM1'])),abs(g[1]-p[1]-float(r['delta_coding_log2TPM1'])))
assert len(samples)==66 and len({r['patient'] for r in samples})==14 and len(deltas)==18 and error<1e-10
summary=[]
for site in ['peritoneal','ovary','liver','lung']:
 rr=[r for r in deltas if r['site']==site];v=[float(r['delta_log2TPM1']) for r in rr];cv=[float(r['delta_coding_log2TPM1']) for r in rr];summary.append({'site':site,'patients':len(rr),'mean':statistics.mean(v),'median':statistics.median(v),'positive':sum(t>0 for t in v),'negative':sum(t<0 for t in v),'coding_mean':statistics.mean(cv),'coding_positive':sum(t>0 for t in cv),'coding_negative':sum(t<0 for t in cv),'leave_one_out_mean_range':[min(statistics.mean(v[:i]+v[i+1:]) for i in range(len(v))),max(statistics.mean(v[:i]+v[i+1:]) for i in range(len(v)))] if len(v)>1 else None})
organ=[]
for donor in sorted({r['patient'] for r in samples}):
 if (donor,'peritoneal') in means and (donor,'ovary') in means:organ.append({'patient':donor,'PM_minus_ovary':means[(donor,'peritoneal')][0]-means[(donor,'ovary')][0],'coding_PM_minus_ovary':means[(donor,'peritoneal')][1]-means[(donor,'ovary')][1]})
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'independent_csv_scalar_recalculation_max_error':error,'samples_verified':len(samples),'patients':14,'site_pairs':18,'route_summary':summary,'within_patient_organ_contrast':organ,'organ_contrast_role':'post-summarydescriptivecontrast;noinferentialP,notindependentvalidation','peritoneal_patient_values':[x for x in deltas if x['site']=='peritoneal']};(O/'results/ROUTE_NUMERICAL_QA_100.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
