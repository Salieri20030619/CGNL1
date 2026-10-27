import pathlib,re,gzip
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/01_metadata')
for f in ['iprox_common.js','iprox_MSV02201.js','iprox_MSV02202.js']:
 t=(r/f).read_text();print('FILE',f)
 for key in ['jsonServer','function project(','findBySubProjectId','findByProjectId']:
  m=t.find(key)
  if m>=0:print(key,t[max(0,m-100):m+500])
for l in gzip.open(r/'GSE314812.soft.gz','rt'):
 if l.startswith(('^SAMPLE','!Sample_title','!Sample_source_name','!Sample_description')):print(l.strip())
