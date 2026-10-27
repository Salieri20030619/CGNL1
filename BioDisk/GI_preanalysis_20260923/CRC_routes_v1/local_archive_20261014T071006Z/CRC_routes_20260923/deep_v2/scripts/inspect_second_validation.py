import pathlib,gzip,re,io,pandas as pd
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');s=gzip.open(O/'metadata/GSE14297.soft.gz','rt',encoding='utf-8').read();p=s.split('!platform_table_begin\n')[1].split('!platform_table_end')[0];df=pd.read_csv(io.StringIO(p),sep='\t',low_memory=False);print('PLATFORM',df.shape,df.columns.tolist());print(df.head(2).to_string(index=False));df.to_csv(O/'metadata/GPL6370_annotation.tsv.gz',sep='\t',index=False)
rows=[]
for c in s.split('^SAMPLE = ')[1:]:
 acc=c.splitlines()[0];title=re.search(r'!Sample_title = (.+)',c).group(1);chars=re.findall(r'!Sample_characteristics_ch1 = (.+)',c);print(acc,title,chars[:1]);rows.append({'sample':acc,'title':title,'characteristics':';'.join(chars)})
pd.DataFrame(rows).to_csv(O/'metadata/GSE14297_samples.tsv',sep='\t',index=False)
c=s.split('^SAMPLE = ')[1];print('PROCESS','\n'.join(x for x in c.splitlines() if x.startswith('!Sample_data_processing')))
with gzip.open(O/'metadata/GSE14297.matrix.gz','rt') as f:
 for line in f:
  if line.startswith('!series_matrix_table_begin'):
   print('MATRIX',next(f)[:250],next(f)[:250]);break
