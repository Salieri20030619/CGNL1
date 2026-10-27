from pathlib import Path
import gzip,csv,json,datetime
import numpy as np
p=Path.cwd();ids=[];values=[];header=None;inside=False
with gzip.open(p/'raw/GSE140973_series_matrix.txt.gz','rt',encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!series_matrix_table_begin':inside=True;continue
  if row[0]=='!series_matrix_table_end':break
  if not inside:continue
  if header is None:header=row;assert header[0]=='ID_REF';continue
  assert len(row)==len(header);ids.append(row[0]);values.append([float(x) if x not in ['','null','NA','NaN','nan','N/A'] else np.nan for x in row[1:]])
a=np.asarray(values,dtype=float);assert len(set(ids))==len(ids);assert len(set(header[1:]))==114;meta=json.loads((p/'metadata/SAMPLE_PAIRING_296b.json').read_text());assert set(header[1:])=={r['accession'] for r in meta['samples']}
finite=np.isfinite(a);assert not np.isinf(a).any();quant=np.nanquantile(a,[0,.01,.25,.5,.75,.99,1],axis=0);complete=a[finite.all(axis=1)];sorted_complete=np.sort(complete,axis=0);maxsort=float(np.max(np.abs(sorted_complete-sorted_complete[:,[0]])))
np.savez_compressed(p/'processed/matrix_input_306.npz',probe_ids=np.asarray(ids),accessions=np.asarray(header[1:]),values=a)
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'probes':len(ids),'samples':len(header)-1,'missing_values':int((~finite).sum()),'negative_values':int((a<0).sum()),'zero_values':int((a==0).sum()),'complete_probes':len(complete),'global_range':[float(np.nanmin(a)),float(np.nanmax(a))],'quantile_column_ranges':{str(q):[float(min(v)),float(max(v))] for q,v in zip([0,.01,.25,.5,.75,.99,1],quant)},'max_sorted_complete_column_difference':maxsort,'target_expression_not_inspected':True,'preprocessing_decision_pending':True}
(p/'metadata/MATRIX_SCALE_QA_306.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
