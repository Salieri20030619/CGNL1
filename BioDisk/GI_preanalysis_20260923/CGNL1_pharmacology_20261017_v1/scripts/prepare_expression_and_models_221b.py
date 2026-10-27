
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import csv,zipfile,io,json,datetime
import numpy as np,pandas as pd
p=Path.cwd();genes=['CGNL1','EPCAM','VIM','MKI67'];measurements={};geneids={}
with zipfile.ZipFile(p/'raw/rnaseq_20191101.zip') as z:
 for member,label in [('rnaseq_read_count_20191101.csv','counts'),('rnaseq_fpkm_20191101.csv','FPKM')]:
  with io.TextIOWrapper(z.open(member),encoding='utf-8') as h:
   rd=csv.reader(h);header=next(rd);names=next(rd);sources=next(rd);ann=next(rd)
   model_ids=header[2:];assert len(set(model_ids))==len(model_ids)
   values={};total=np.zeros(len(model_ids));gene_rows=0;common_gene_ids=[]
   for row in rd:
    assert len(row)==len(header)
    if label=='counts':
     if '' not in row[2:]:
      v=np.array(row[2:],dtype=float);assert np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all();total+=v;gene_rows+=1;common_gene_ids.append(row[0])
    if row[1] in genes:
     assert row[1] not in values;values[row[1]]=np.array(row[2:],dtype=float);geneids[row[1]]=row[0]
   assert set(values)==set(genes)
   if label=='counts':
    assert (total>0).all() and gene_rows==34985
    (p/'metadata/common_denominator_gene_ids_221b.json').write_text(json.dumps(common_gene_ids),encoding='utf-8')
    frame=pd.DataFrame({'model_id':model_ids,'RNA_model_name':names[2:],'RNA_dataset_source':sources[2:],'common_gene_read_count_sum':total})
    for g in genes:frame[g+'_count']=values[g];frame[g+'_log2CPM1']=np.log2(1+values[g]/total*1e6)
   else:
    assert model_ids==frame.model_id.tolist()
    for g in genes:frame[g+'_FPKM']=values[g];frame[g+'_log2FPKM1']=np.log2(1+values[g])
assert geneids['CGNL1']=='SIDG04503'
models=pd.read_csv(p/'raw/model_list_20230608.csv',dtype=str,keep_default_na=False)
assert models.model_id.is_unique
cols=['model_id','patient_id','sample_id','model_name','tissue','cancer_type','cancer_type_detail','model_type','species','COSMIC_ID','CCLE_ID']
joined=frame.merge(models[cols],on='model_id',how='left',validate='one_to_one',indicator=True)
assert (joined['_merge']=='both').all()
types={'Pancreatic Carcinoma':'Pancreatic','Colorectal Carcinoma':'Colorectal','Hepatocellular Carcinoma':'Hepatocellular','Biliary Tract Carcinoma':'Biliary','Gastric Carcinoma':'Gastric','Esophageal Squamous Cell Carcinoma':'Esophageal_squamous','Esophageal Carcinoma':'Esophageal_unspecified'}
tissues=['Pancreas','Large Intestine','Small Intestine','Liver','Biliary Tract','Stomach','Esophagus']
joined['eligible']=joined.tissue.isin(tissues)&joined.cancer_type.isin(types)&joined.model_type.eq('Cell Line')&joined.species.eq('Homo Sapiens')&joined.patient_id.ne('')
joined['lineage']=joined.cancer_type.map(types).fillna('outside_scope')
mask=joined.cancer_type.eq('Esophageal Carcinoma')&joined.cancer_type_detail.str.contains('adenocarcinoma',case=False)
joined.loc[mask,'lineage']='Esophageal_adenocarcinoma'
eligible=joined[joined.eligible].sort_values('model_id').copy()
selected=eligible.drop_duplicates('patient_id',keep='first').copy()
excluded=eligible[~eligible.model_id.isin(selected.model_id)].copy()
assert selected.patient_id.is_unique and selected.model_id.is_unique
frame.to_csv(p/'processed/expression_target_covariates_all_models_221.tsv',sep=chr(9),index=False,encoding='utf-8')
joined.to_csv(p/'processed/model_eligibility_221.tsv',sep=chr(9),index=False,encoding='utf-8')
selected.to_csv(p/'processed/selected_patient_models_221.tsv',sep=chr(9),index=False,encoding='utf-8')
excluded.to_csv(p/'processed/same_patient_excluded_models_221.tsv',sep=chr(9),index=False,encoding='utf-8')
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'denominator':'34985 source-common measured genes','RNA_models':len(frame),'eligible_GI_models':len(eligible),'selected_unique_patients':len(selected),'same_patient_models_excluded':len(excluded),'lineages':selected.lineage.value_counts().to_dict(),'RNA_sources':selected.RNA_dataset_source.value_counts().to_dict(),'source_gene_ids':geneids,'target_nonzero_count':int((selected.CGNL1_count>0).sum()),'target_integer_FPKM_zero':int((selected.CGNL1_FPKM==0).sum()),'target_log2CPM1_range':[float(selected.CGNL1_log2CPM1.min()),float(selected.CGNL1_log2CPM1.max())],'response_associations_not_inspected':True}
(p/'metadata/MODEL_MAPPING_QUALIFICATION_221.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
