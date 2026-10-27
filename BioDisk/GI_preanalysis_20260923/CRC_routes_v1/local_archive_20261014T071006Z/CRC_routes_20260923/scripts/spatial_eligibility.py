import pathlib,pandas as pd,json,gzip,subprocess
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');m=pd.read_csv(O/'metadata/spatial_author_metadata.tsv.gz',sep='\t');m['sample']=m.Barcode.str.split('~').str[0];z=m[m.Layer1=='Primary Colon'];print(pd.crosstab(z['sample'],z.Layer3).to_string());z.groupby(['sample','Layer3']).size().rename('spots').reset_index().to_csv(O/'results/spatial_region_coverage.tsv',sep='\t',index=False)
j=json.loads((O/'metadata/spatial_table2.json').read_text(encoding='utf-8'));rows=[];patient=None;stype=None
for r in j['Tab 1'][3:]:
 if isinstance(r[0],int):patient=r[0];stype=None
 if r[2]:stype=r[2]
 if r[1] and patient:rows.append({'patient':patient,'sample':r[1].replace('IX','IX_') if r[1].startswith('IX') else r[1],'sample_type':stype,'visium':r[3]=='✓'})
pd.DataFrame(rows).to_csv(O/'metadata/spatial_patient_map.tsv',sep='\t',index=False)
q=z[['sample']].drop_duplicates().merge(pd.DataFrame(rows),on='sample');print(q.to_string(index=False))
