
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,itertools
import numpy as np,pandas as pd
from scipy import stats
p=Path.cwd();spec=json.loads((p/'STATISTICAL_SPECIFICATION_227.json').read_text());assert json.loads((p/'metadata/GSE210993_DOWNLOAD_COMPLETE_227.json').read_text())['success']
samples=pd.DataFrame(json.loads((p/'metadata/GSE210993_ANALYSIS_SAMPLES_227.json').read_text()))
x=pd.read_csv(p/'raw/GSE210993_IM.expression_TPM_matrix.103.txt.gz',sep=chr(9))
names=samples.sample_name.tolist();assert len(names)==103 and len(set(names))==103 and all(n in x.columns for n in names)
a=x[names].to_numpy(float);assert np.isfinite(a).all() and (a>=0).all();sums=a.sum(axis=0);assert np.max(abs(sums-1e6))<5000
target=x[x.Gene_ID.astype(str).str.split('.').str[0]=='ENSG00000128849'];assert len(target)==1 and str(target.Gene_Name.iloc[0]).strip('"')=='CGNL1'
samples['CGNL1_TPM']=[float(target[n].iloc[0]) for n in names];samples['CGNL1_log2TPM1']=np.log2(1+samples.CGNL1_TPM)
context=[]
for g in spec['context_genes']:
 row=x[x.Gene_Name==g];assert len(row)==1,(g,len(row))
 for n in names:context.append({'sample_name':n,'gene':g,'TPM':float(row[n].iloc[0]),'log2TPM1':float(np.log2(1+row[n].iloc[0]))})
summary=[];deltas=[]
def compare(label,d,role,fdr_family=None):
 d=np.asarray(d,dtype=float);rec={'contrast':label,'role':role,'family':fdr_family,'independent_patients':len(d),'mean_difference_log2TPM1':float(d.mean()) if len(d) else None,'median_difference_log2TPM1':float(np.median(d)) if len(d) else None,'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P':None,'FDR':None}
 if len(d)>=3 and (d!=0).any():
  nz=d[d!=0];ranks=stats.rankdata(abs(nz));obs=abs(np.dot(np.sign(nz),ranks))
  rec['P']=float(np.mean([abs(np.dot(sign,ranks))>=obs-1e-12 for sign in itertools.product([-1,1],repeat=len(nz))]))
 summary.append(rec)
for sensitivity in [False,True]:
 s=samples if sensitivity else samples[~samples.reclassified]
 s=s[s.type.isin(['Human IM organoid','Human normal gastric organoid'])]
 w=s.groupby(['patient_id','type']).CGNL1_log2TPM1.mean().unstack().dropna()
 diff=w['Human IM organoid']-w['Human normal gastric organoid']
 if not sensitivity:assert sorted(diff.index.tolist())==spec['paired_patient_ids']
 label='IM_organoid_minus_normal_gastric'+('_author_reclassification_sensitivity' if sensitivity else '_primary')
 compare(label,diff,'sensitivity' if sensitivity else 'primary')
 for patient,value in diff.items():deltas.append({'contrast':label,'patient_id':patient,'difference_log2TPM1':float(value)})
idx=samples.set_index('sample_name')
longpairs=[('GX001-IM3O_P30','GX001-IM3O'),('GX059-IMOL','GX059-IMO-P15'),('GX072-IM3OL','GX072-IM3O'),('GX072-IMOL','GX072-IMO'),('GX076-IMOL','GX076-IMO')]
pair_records=[]
for late,early in longpairs:
 assert idx.loc[late,'patient_id']==idx.loc[early,'patient_id']
 pair_records.append({'pair_type':'long_term','patient_id':idx.loc[late,'patient_id'],'case_sample':late,'reference_sample':early,'difference_log2TPM1':float(idx.loc[late,'CGNL1_log2TPM1']-idx.loc[early,'CGNL1_log2TPM1'])})
long=pd.DataFrame(pair_records).groupby('patient_id').difference_log2TPM1.mean()
compare('long_term_minus_baseline_IM',long,'secondary','secondary_target_contrasts')
for patient,value in long.items():deltas.append({'contrast':'long_term_minus_baseline_IM','patient_id':patient,'difference_log2TPM1':float(value)})
unmatched=[]
for _,r in samples[(samples.type=='Human IM organoid')&~samples.reclassified].iterrows():
 name=r.sample_name
 base=name[:-4] if name.endswith('-P15') else name
 tissue=base[:-1]+'F' if base.endswith('O') else None
 if tissue not in idx.index or idx.loc[tissue,'type']!='Human IM frozen tissue':unmatched.append(name);continue
 assert idx.loc[tissue,'patient_id']==r.patient_id
 pair_records.append({'pair_type':'organoid_minus_tissue','patient_id':r.patient_id,'case_sample':name,'reference_sample':tissue,'difference_log2TPM1':float(r.CGNL1_log2TPM1-idx.loc[tissue,'CGNL1_log2TPM1'])})
pr=pd.DataFrame(pair_records);tissue=pr[pr.pair_type=='organoid_minus_tissue'].groupby('patient_id').difference_log2TPM1.mean()
compare('IM_organoid_minus_matched_IM_tissue',tissue,'secondary','secondary_target_contrasts')
for patient,value in tissue.items():deltas.append({'contrast':'IM_organoid_minus_matched_IM_tissue','patient_id':patient,'difference_log2TPM1':float(value)})
res=pd.DataFrame(summary);mask=res.role=='secondary';valid=mask&res.P.notna();res.loc[valid,'FDR']=stats.false_discovery_control(res.loc[valid,'P'].to_numpy())
out=p/'results/GSE210993_bulk_230';out.mkdir(exist_ok=True)
samples.to_csv(out/'target_samples.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(context).to_csv(out/'context_gene_samples.tsv',sep=chr(9),index=False,encoding='utf-8')
pr.to_csv(out/'matched_specimens.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(deltas).to_csv(out/'patient_deltas.tsv',sep=chr(9),index=False,encoding='utf-8')
res.to_csv(out/'target_contrasts.tsv',sep=chr(9),index=False,encoding='utf-8')
descriptive=samples.groupby('type').agg(profiles=('sample_name','size'),patients=('patient_id','nunique'),median_TPM=('CGNL1_TPM','median'),min_TPM=('CGNL1_TPM','min'),max_TPM=('CGNL1_TPM','max'))
descriptive.to_csv(out/'group_description.tsv',sep=chr(9),encoding='utf-8')
done={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'profiles':103,'TPM_sum_range':[float(sums.min()),float(sums.max())],'target_gene_id':'ENSG00000128849','contrasts':res.replace({np.nan:None}).to_dict('records'),'unmatched_organoid_tissue_samples':unmatched,'independent_numerical_QA_pending':True,'mechanism_established':False}
(out/'ANALYSIS_COMPLETE_230.json').write_text(json.dumps(done,indent=2),encoding='utf-8')
print(json.dumps(done,indent=2));print(descriptive.to_string())
