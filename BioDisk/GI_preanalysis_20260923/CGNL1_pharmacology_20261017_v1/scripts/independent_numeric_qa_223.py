
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import csv,zipfile,io,json,math,datetime
import numpy as np,pandas as pd
from scipy.stats import t
p=Path.cwd();models=pd.read_csv(p/'processed/selected_patient_models_221.tsv',sep=chr(9));assert models.patient_id.is_unique and models.model_id.is_unique
geneids=set(json.loads((p/'metadata/common_denominator_gene_ids_221b.json').read_text()))
assert len(geneids)==34985
original=models.set_index('model_id');names=models.model_id.tolist();genes=['CGNL1','EPCAM','VIM','MKI67'];totals=[0]*len(names);target={}
with zipfile.ZipFile(p/'raw/rnaseq_20191101.zip') as z:
 with io.TextIOWrapper(z.open('rnaseq_read_count_20191101.csv'),encoding='utf-8') as h:
  rd=csv.reader(h);header=next(rd);next(rd);next(rd);next(rd);index=[header.index(n) for n in names];seen=0
  for row in rd:
   if row[0] in geneids:
    vals=[int(float(row[k])) for k in index]
    totals=[a+b for a,b in zip(totals,vals)];seen+=1
   if row[1] in genes:target[row[1]]=[int(float(row[k])) for k in index]
assert seen==34985 and set(target)==set(genes)
maxnorm=0.
for i,name in enumerate(names):
 assert totals[i]==original.loc[name,'common_gene_read_count_sum']
 for g in genes:
  assert target[g][i]==original.loc[name,g+'_count']
  error=abs(math.log2(1+target[g][i]/totals[i]*1e6)-original.loc[name,g+'_log2CPM1']);assert error<1e-12;maxnorm=max(maxnorm,error)
results=pd.read_csv(p/'results/association_models_222.tsv',sep=chr(9));checks=[]
for _,r in results[results.status=='evaluated'].iterrows():
 f=p/'processed/model_fits_222'/(r.screen+'_'+r.drug+'_'+r.role+'.tsv')
 d=pd.read_csv(f,sep=chr(9));assert d.patient_id.is_unique and d.model_id.is_unique
 cols=[c for c in d.columns if c not in ['response','patient_id','model_id']]
 X=d[cols].to_numpy(float);y=d.response.to_numpy(float);assert len(d)==r.n_patient_models
 inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;res=y-X@beta;h=np.sum((X@inv)*X,axis=1)
 omega=(res/(1-h))**2;cov=inv@((X.T*omega)@X)@inv
 j=cols.index('CGNL1');se=math.sqrt(cov[j,j]);df=len(y)-X.shape[1];pval=2*t.sf(abs(beta[j]/se),df);low=beta[j]-t.ppf(.975,df)*se
 errors={'beta':abs(beta[j]-r.beta_CGNL1),'SE':abs(se-r.SE_HC3),'P':abs(pval-r.P),'CI_low':abs(low-r.CI_low)}
 assert max(errors.values())<1e-8,(f.name,errors)
 checks.append({'model':f.stem,'n':len(d),'max_error':max(errors.values())})
maxfdr=0.
for keys,group in results.groupby(['screen','family','role','response','expression_measure']):
 v=group.dropna(subset=['P']).sort_values('P');n=len(v)
 if not n:continue
 scaled=v.P.to_numpy()*n/np.arange(1,n+1);q=np.minimum.accumulate(scaled[::-1])[::-1];q=np.minimum(q,1.)
 error=float(np.max(np.abs(q-v.FDR.to_numpy())));assert error<1e-12;maxfdr=max(maxfdr,error)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'selected_unique_patient_models':len(models),'normalization_source_common_genes':seen,'raw_integer_count_and_normalization_max_error':maxnorm,'regression_models_checked':len(checks),'manual_HC3_matrix_max_error':max(c['max_error'] for c in checks),'manual_BH_max_error':maxfdr,'checks':checks,'scope':'all172 selected models count normalization; all evaluated regression designs, target coefficients, HC3 SE, P and lowerCI; all reported FDR families','independent_clinical_validation':False}
(p/'results/INDEPENDENT_NUMERICAL_QA_223.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
