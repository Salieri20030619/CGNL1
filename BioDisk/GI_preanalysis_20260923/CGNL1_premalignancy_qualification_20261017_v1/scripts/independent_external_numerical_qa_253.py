from pathlib import Path
import json,csv,gzip,math,datetime
import numpy as np
from scipy.stats import t
p=Path.cwd();out=p/'results/external_252';a=json.loads((out/'ANALYSIS_COMPLETE_252.json').read_text());probes=a['diagnostic_probes'];primary_probes=a['eligible_probes'];raw={};header=None
with gzip.open(p/'raw/external_250/GSE78523_series_matrix.txt.gz','rt',encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row or row[0].startswith('!'):continue
  if header is None:header=row[1:];continue
  if row[0] in probes:
   assert row[0] not in raw;raw[row[0]]={g:float(v) for g,v in zip(header,row[1:])}
with (out/'sample_source_values.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==45 and len({r['accession'] for r in rows})==45
for r in rows:
 for probe in probes:assert abs(float(r[probe])-raw[probe][r['accession']])<1e-12
 assert abs(float(r['CGNL1_mean_log2'])-math.fsum(raw[probe][r['accession']] for probe in primary_probes)/len(primary_probes))<1e-12
metadata=json.loads((p/'processed/external_251/sample_metadata_rows.json').read_text())
accessions=next(r[1:] for r in metadata if r[0]=='!Sample_geo_accession');characteristics=[r[1:] for r in metadata if r[0]=='!Sample_characteristics_ch1'];source={g:{} for g in accessions}
for field in characteristics:
 for g,value in zip(accessions,field):
  key,val=value.split(': ',1);source[g][key.lower()]=val
for r in rows:
 for key in ['status','sex','localization']:assert r[key]==source[r['accession']][key]
 assert float(r['age'])==float(source[r['accession']]['age'])
maxerr=0.
for result in a['all_results']:
 with (out/(result['contrast']+'_model_values.tsv')).open(encoding='utf-8',newline='') as f:mr=list(csv.DictReader(f,delimiter='\t'))
 label=result['contrast'];expected={g for g,v in source.items() if (v['status'].startswith('IIM') if label=='secondary_IIM_progression' else v['status'].startswith('CIM') if label=='secondary_CIM_progression' else not v['status'].endswith('-GC') if label=='secondary_IM_vs_healthy' else v['status']!='Healthy Control')}
 assert {r['accession'] for r in mr}==expected
 for row in mr:
  v=source[row['accession']]
  for key in ['progression','IIM','IM','man']:
   if key in row:
    expected_value={'progression':v['status'].endswith('-GC'),'IIM':v['status'].startswith('IIM'),'IM':v['status']!='Healthy Control','man':v['sex']=='Man'}[key];assert float(row[key])==float(expected_value)
  if 'age' in row:assert abs(float(row['age'])-(float(v['age'])-math.fsum(float(source[g]['age']) for g in expected)/len(expected)))<1e-10
 cols=result['model_columns'].split(',');X=np.array([[float(r[k]) for k in cols] for r in mr]);y=np.array([float(r[result['outcome']]) for r in mr]);n,k=X.shape;inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;res=y-X@beta;h=np.array([x@inv@x for x in X]);assert max(h)<1
 meat=np.zeros((k,k))
 for xi,ei,hi in zip(X,res,h):meat+=np.outer(xi,xi)*(ei/(1-hi))**2
 cov=inv@meat@inv;j=cols.index(result['target']);se=math.sqrt(cov[j,j]);pval=2*t.sf(abs(beta[j]/se),n-k);crit=t.ppf(.975,n-k)
 vals={'effect':beta[j],'CI_low':beta[j]-crit*se,'CI_high':beta[j]+crit*se,'P':pval}
 for key,val in vals.items():
  err=abs(float(val)-result[key]);maxerr=max(maxerr,err);assert err<1e-8,(result['contrast'],key,err)
for family in ['secondary','probe_sensitivity']:
 rr=[r for r in a['all_results'] if r['role']==family];order=sorted(range(len(rr)),key=lambda i:rr[i]['P']);adjusted=[None]*len(rr);running=1.
 for rank in range(len(rr),0,-1):
  i=order[rank-1];running=min(running,rr[i]['P']*len(rr)/rank);adjusted[i]=running
 for r,v in zip(rr,adjusted):assert abs(r['FDR']-v)<1e-12
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'source_samples_checked':45,'diagnostic_probes_checked':len(probes),'primary_probes_checked':len(primary_probes),'models_checked':len(a['all_results']),'max_fit_absolute_error':maxerr,'methods':'Raw CSV source values and scalar probe means; independent normal-equation estimates, explicit HC3 sandwich and t intervals; separate manual BH families','whole_project_complete':False}
(out/'INDEPENDENT_NUMERICAL_QA_253.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
