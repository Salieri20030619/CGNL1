from pathlib import Path
import csv,io,zipfile,json,math,datetime,collections
import numpy as np
p=Path.cwd();out=p/'results/organoid_264';out.mkdir(exist_ok=False)
assert json.loads((p/'metadata/DOWNLOAD_COMPLETE_263.json').read_text())['all_pass']
with (p/'raw/model_list_20261013.csv').open(encoding='utf-8') as f:models=list(csv.DictReader(f))
lookup={r['model_id']:r for r in models};assert len(lookup)==len(models)
data={};header_ref=None;genes_ref=None;target_ids={'CGNL1':'ENSG00000128849','CGN':'ENSG00000143375','CDH1':'ENSG00000039068'};target_rows={};target_annotations={}
with zipfile.ZipFile(p/'raw/Organoids_fitness_scores_20260723.zip') as z:
 for name in z.namelist():
  key='binary' if 'binary_matrix' in name else 'scaled_BF' if 'scaled_bayesian' in name else 'BF' if 'bayesian_factors' in name else 'fold_change'
  with z.open(name) as f:
   reader=csv.reader(io.TextIOWrapper(f,encoding='utf-8'),delimiter='\t');headers=[next(reader) for _ in range(5)];assert [r[0] for r in headers]==['model_name','model_id','source','qc_pass','gene_id'];genes=[];values=[];targets={}
   for row in reader:
    assert len(row)==len(headers[0]);genes.append(tuple(row[:3]));v=[float(x) if x not in ['', 'NA', 'NaN','nan'] else math.nan for x in row[3:]];values.append(v)
    if row[1] in target_ids:
     assert row[2]==target_ids[row[1]] and row[1] not in targets;targets[row[1]]=v;target_annotations[row[1]]=row[:3]
  if header_ref is None:header_ref=headers;genes_ref=genes
  else:assert headers==header_ref and genes==genes_ref
  assert set(targets)==set(target_ids);data[key]=np.array(values);target_rows[key]=targets
assert len(set(genes_ref))==len(genes_ref)
mask=np.isfinite(data['scaled_BF']) & np.isfinite(data['binary']);bv=data['binary'][mask];assert set(np.unique(bv))<= {0.,1.}
agree_gt=float(np.mean((data['scaled_BF'][mask]>0)==bv));agree_lt=float(np.mean((data['scaled_BF'][mask]<0)==bv));assert max(agree_gt,agree_lt)==1.,(agree_gt,agree_lt)
records=[]
for i,mid in enumerate(header_ref[1][3:]):
 assert mid in lookup;meta=lookup[mid];r={k:meta[k] for k in ['model_id','model_name','patient_id','sample_id','tissue','cancer_type','tissue_status','model_type','parent_id']};r['column_index']=i;r['source']=header_ref[2][i+3];r['qc_pass']=header_ref[3][i+3]
 for measure in data:
  for gene in target_ids:r[gene+'_'+measure]=target_rows[measure][gene][i]
 records.append(r)
with (out/'target_model_values.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
qc=[r for r in records if r['qc_pass'].upper()=='TRUE'];eligible=[r for r in qc if math.isfinite(r['CGNL1_binary'])];dependent=[r for r in eligible if r['CGNL1_binary']==1]
by=[]
for tissue in sorted({r['tissue'] for r in eligible}):
 rr=[r for r in eligible if r['tissue']==tissue];dd=[r for r in rr if r['CGNL1_binary']==1];by.append({'tissue':tissue,'models':len(rr),'unique_patients':len({r['patient_id'] for r in rr if r['patient_id']}),'CGNL1_dependent_models':len(dd)})
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matrix_columns':len(records),'unique_models':len({r['model_id'] for r in records}),'genes':len(genes_ref),'qc_pass_columns':len(qc),'CGNL1_evaluable_columns':len(eligible),'unique_patients':len({r['patient_id'] for r in eligible if r['patient_id']}),'missing_patient_ids':sum(not r['patient_id'] for r in eligible),'target_annotations':target_annotations,'binary_matches_positive_scaled_BF':agree_gt,'binary_matches_negative_scaled_BF':agree_lt,'all_matrix_binary_checks':int(mask.sum()),'CGNL1_dependent_models':dependent,'by_tissue':by,'source_call_definition':'Source BAGEL-derived binary screen call, not newly computed biological P; scaled-BF file sign verified against binary matrix.','association_tests_run':False,'whole_project_complete':False}
(out/'QUALIFICATION_COMPLETE_264.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
