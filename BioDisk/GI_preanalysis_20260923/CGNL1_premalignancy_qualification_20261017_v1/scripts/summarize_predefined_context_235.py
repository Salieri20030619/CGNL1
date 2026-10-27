
from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
p=Path.cwd();out=p/'results/GSE210993_bulk_230'
assert json.loads((out/'INDEPENDENT_NUMERICAL_QA_232.json').read_text())['pass']
meta=pd.DataFrame(json.loads((p/'metadata/GSE210993_ANALYSIS_SAMPLES_227.json').read_text()))
context=pd.read_csv(out/'context_gene_samples.tsv',sep=chr(9));target=pd.read_csv(out/'target_samples.tsv',sep=chr(9))
t=target[['sample_name','CGNL1_TPM','CGNL1_log2TPM1']].rename(columns={'CGNL1_TPM':'TPM','CGNL1_log2TPM1':'log2TPM1'});t['gene']='CGNL1';t['measurement_status']='measured'
allgenes=pd.concat([context,t],ignore_index=True).merge(meta[['sample_name','patient_id','type','reclassified']],on='sample_name',validate='many_to_one')
patients=json.loads((p/'STATISTICAL_SPECIFICATION_227.json').read_text())['paired_patient_ids'];s=allgenes[~allgenes.reclassified&allgenes.patient_id.isin(patients)]
rows=[]
for gene,g in s.groupby('gene'):
 g=g[g.type.isin(['Human IM organoid','Human normal gastric organoid'])]
 w=g.groupby(['patient_id','type']).log2TPM1.mean().unstack().dropna()
 d=w['Human IM organoid']-w['Human normal gastric organoid'] if len(w) else pd.Series(dtype=float)
 rows.append({'gene':gene,'independent_pairs':len(d),'mean_delta_log2TPM1':float(d.mean()) if len(d) else np.nan,'median_delta_log2TPM1':float(d.median()) if len(d) else np.nan,'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P_computed':False,'role':'predefined_context_descriptive'})
r=pd.DataFrame(rows);r.to_csv(out/'paired_context_description_235.tsv',sep=chr(9),index=False,encoding='utf-8')
patientlevels=allgenes.groupby(['gene','patient_id','type'],as_index=False)[['TPM','log2TPM1']].mean()
patientlevels.to_csv(out/'context_patient_level_values_235.tsv',sep=chr(9),index=False,encoding='utf-8')
summary=patientlevels.groupby(['gene','type']).agg(patients=('patient_id','size'),median_patient_mean_TPM=('TPM','median')).reset_index()
summary.to_csv(out/'context_tissue_class_description_235.tsv',sep=chr(9),index=False,encoding='utf-8')
print(r[r.gene.isin(['CGNL1','CGN','CDH1','TJP1','CAMSAP3','EPCAM','MKI67','CDX2','VIL1','MUC2','MUC5AC','GKN1'])].to_string(index=False))
print(summary[summary.gene.isin(['CGNL1','CGN'])].to_string(index=False))
