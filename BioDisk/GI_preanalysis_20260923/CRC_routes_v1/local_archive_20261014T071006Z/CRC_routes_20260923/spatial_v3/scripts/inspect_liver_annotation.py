import pathlib,pandas as pd,json,gzip,re,sys
sys.stdout.reconfigure(encoding='utf-8')
B=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');O=B/'spatial_v3';m=pd.read_csv(B/'metadata/spatial_author_metadata.tsv.gz',sep='\t');m['sample']=m.Barcode.str.split('~').str[0];pm=pd.read_csv(B/'metadata/spatial_patient_map.tsv',sep='\t');p=pm[pm.visium & pm.sample_type.str.contains('Liver')];a=m[m['sample'].isin(p['sample'])];print('COLUMNS',list(m.columns));print('LAYER1',a.Layer1.value_counts().to_dict());print('LAYER2',a.Layer2.value_counts().to_dict());print('LAYER3',a.Layer3.value_counts().to_dict());print(pd.crosstab(a['sample'],a.Layer3).to_string());print('EXAMPLE',a.head(2).to_dict(orient='records'))
p.merge(a.groupby('sample').size().rename('annotated_spots'),on='sample',how='left').to_csv(O/'results/liver_annotation_eligibility.tsv',sep='\t',index=False)
s=gzip.open(B/'metadata/GSE294385.soft.gz','rt',encoding='utf-8').read();jobs=[]
for c in s.split('^SAMPLE = ')[1:]:
 title=re.search(r'!Sample_title = (.+)',c).group(1).strip()
 if title in set(p['sample']):
  u=re.search(r'!Sample_supplementary_file_\d+ = (ftp.+)',c).group(1).strip().replace('ftp://','https://');jobs.append({'sample':title,'url':u,'patient':int(p.set_index('sample').loc[title,'patient']),'sample_type':p.set_index('sample').loc[title,'sample_type']})
(O/'metadata/LIVER_DOWNLOAD_PLAN.json').write_text(json.dumps(jobs,indent=2),encoding='utf-8');print('DOWNLOADS',len(jobs));print(jobs[0])
