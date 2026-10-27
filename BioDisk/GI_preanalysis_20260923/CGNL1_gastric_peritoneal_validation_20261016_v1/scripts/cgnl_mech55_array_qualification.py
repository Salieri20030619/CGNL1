from pathlib import Path
import json,io
import pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';s=(O/'metadata/GSE53276_samples.soft').read_bytes().decode('utf-8');rows=[]
for b in s.split('^SAMPLE = ')[1:]:
 ls=b.splitlines();beg=[i for i,l in enumerate(ls) if l.strip()=='!sample_table_begin'];end=[i for i,l in enumerate(ls) if l.strip()=='!sample_table_end'];r={'accession':ls[0],'begin':beg,'end':end}
 if beg and end:
  tab=pd.read_csv(io.StringIO('\n'.join(ls[beg[0]+1:end[0]])),sep='\t');r.update(shape=list(tab.shape),columns=list(tab.columns),first_row=tab.iloc[0].to_dict())
 rows.append(r)
print(json.dumps(rows,indent=2));p=O/'metadata/GPL570_full.soft'
if p.exists():
 ls=p.read_bytes().decode('utf-8').splitlines();beg=[i for i,l in enumerate(ls) if l.strip()=='!platform_table_begin'];end=[i for i,l in enumerate(ls) if l.strip()=='!platform_table_end'];print('GPL markers',beg,end)
 if beg and end:
  tab=pd.read_csv(io.StringIO('\n'.join(ls[beg[0]+1:end[0]])),sep='\t',dtype=str,keep_default_na=False);print('GPL_COLUMNS',list(tab.columns));hit=tab[tab.apply(lambda col:col.str.contains('CGNL1',regex=False)).any(axis=1)];hit.to_csv(O/'metadata/CGNL1_GPL570_probe_annotation.tsv',sep='\t',index=False,encoding='utf-8');print('TARGET_ANNOTATION',hit.to_string(index=False)[:4500])
(O/'metadata/ARRAY_QUALIFICATION_55.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
