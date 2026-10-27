from pathlib import Path
import json,csv,io,zipfile,math,datetime,collections
import numpy as np
from scipy.stats import t
p=Path.cwd();out=p/'results/context_280';a=json.loads((out/'ANALYSIS_COMPLETE_270.json').read_text());symbols={'ENSG00000128849':'CGNL1','ENSG00000143375':'CGN','ENSG00000039068':'CDH1'}
with (p/'raw/model_list_20261013.csv').open(encoding='utf-8') as f:meta={r['model_id']:r for r in csv.DictReader(f)}
with zipfile.ZipFile(p/'raw/Organoids_fitness_scores_20260723.zip') as z:
 name=next(n for n in z.namelist() if 'fold_change_values' in n)
 with z.open(name) as f:
  reader=csv.reader(io.TextIOWrapper(f,encoding='utf-8'),delimiter='\t');headers=[next(reader) for _ in range(5)];mid=headers[1][3:];qc=dict(zip(mid,headers[3][3:]));row=next(r for r in reader if r[2]=='ENSG00000128849');effects=dict(zip(mid,map(float,row[3:])))
patient_models=collections.defaultdict(list);raw_values={}
for m in mid:
 mm=meta[m];dd=json.loads((p/'raw/target_rnaseq_resolved_280'/(m+'.json')).read_text());g={};dataset=set();sources=set()
 for rr in dd['data']:
  v=rr['attributes'];symbol=symbols[v['ensembl_gene_id']];assert symbol not in g;g[symbol]=float(v['tpm']) if v['tpm'] is not None else math.nan;dataset.add(v['dataset_id']);sources.add(v['data_source'])
 raw_values[m]=g
 good=mm['tissue'] in ['Esophagus','Large Intestine','Pancreas','Stomach'] and mm['tissue_status']=='Tumour' and qc[m].upper()=='TRUE' and mm['patient_id'] and set(g)==set(symbols.values()) and all(math.isfinite(v) and v>=0 for v in g.values()) and dataset=={22} and sources=={'Sanger'}
 if good:patient_models[mm['patient_id']].append(m)
with (out/'patient_model_values.tsv').open(encoding='utf-8') as f:pr=list(csv.DictReader(f,delimiter='\t'))
assert {r['patient_id'] for r in pr}==set(patient_models);by_patient={r['patient_id']:r for r in pr}
for r in pr:
 ms=patient_models[r['patient_id']];assert set(r['model_ids'].split(';'))==set(ms);assert int(r['n_models'])==len(ms);assert {meta[m]['tissue'] for m in ms}=={r['tissue']}
 assert abs(float(r['CGNL1_KO_effect'])-math.fsum(effects[m] for m in ms)/len(ms))<1e-12
 for gene in symbols.values():assert abs(float(r[gene+'_log2TPM1'])-math.fsum(raw_values[m][gene] for m in ms)/len(ms))<1e-12
max_error=0.
for res in a['results']:
 gene=res['target'];target=gene+'_log2TPM1'
 with (out/(gene+'_design.tsv')).open(encoding='utf-8') as f:rr=list(csv.DictReader(f,delimiter='\t'))
 assert {r['patient_id'] for r in rr}==set(patient_models);cols=res['columns']
 for row in rr:
  ref=by_patient[row['patient_id']]
  for col in cols:
   expected=1. if col=='intercept' else float(ref['tissue']==col[len('tissue_'):]) if col.startswith('tissue_') else float(ref[col]);assert abs(float(row[col])-expected)<1e-12
  assert float(row['CGNL1_KO_effect'])==float(ref['CGNL1_KO_effect'])
 X=np.array([[float(r[k]) for k in cols] for r in rr]);y=np.array([float(r['CGNL1_KO_effect']) for r in rr]);n,k=X.shape;inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;residual=y-X@beta;h=np.array([v@inv@v for v in X]);assert max(h)<1-1e-10
 meat=np.zeros((k,k))
 for xi,e,hi in zip(X,residual,h):meat+=np.outer(xi,xi)*(e/(1-hi))**2
 cov=inv@meat@inv;j=cols.index(target);se=math.sqrt(cov[j,j]);pv=2*t.sf(abs(beta[j]/se),n-k);cut=t.ppf(.975,n-k)
 for key,value in {'effect':beta[j],'CI_low':beta[j]-cut*se,'CI_high':beta[j]+cut*se,'P':pv}.items():
  err=abs(value-res[key]);max_error=max(max_error,err);assert err<1e-8,(gene,key,err)
 loo=[]
 for i in range(n):
  xx=np.delete(X,i,axis=0);yy=np.delete(y,i);loo.append(float(np.linalg.lstsq(xx,yy,rcond=None)[0][j]))
 assert abs(min(loo)-res['LOO_effect_min'])<1e-8 and abs(max(loo)-res['LOO_effect_max'])<1e-8
rr=a['results'];order=sorted(range(len(rr)),key=lambda i:rr[i]['P']);running=1.
for rank in range(len(rr),0,-1):
 r=rr[order[rank-1]];running=min(running,r['P']*len(rr)/rank);assert abs(r['FDR']-running)<1e-12
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'patients_checked':len(pr),'included_models_checked':sum(len(v) for v in patient_models.values()),'models_checked':2,'max_HC3_absolute_error':max_error,'method':'Independent raw target ZIP row, per-model RNA JSON and clinical metadata; scalar patient aggregation; explicit HC3 normal equations/t inference and manual two-test BH; independent least-squares leave-one-patient-out effects','whole_project_complete':False};(out/'INDEPENDENT_QA_271.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
