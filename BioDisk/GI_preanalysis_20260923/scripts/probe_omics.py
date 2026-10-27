import pathlib,h5py,pandas as pd,gzip
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
p=r/'02_raw'/'PXD064711_protein.raw.intensity.txt'
d=pd.read_csv(p,sep='\t',nrows=3,encoding='utf-16');print('PROTEIN_COLUMNS',d.columns.tolist());print(d.iloc[:2,:15].to_string(index=False))
with h5py.File(r/'02_raw'/'GSE294559_TP001.h5') as h:
 print('H5',list(h),list(h['matrix']),h['matrix']['shape'][:]);print('FEATURE',list(h['matrix']['features']));print('TYPES',pd.Series(h['matrix']['features']['feature_type'][:]).value_counts().to_dict())
for p in (r/'01_metadata').glob('GSE*.soft.gz'):
 if p.name.startswith(('GSE163171','GSE187650','GSE163035','GSE187512')):
  print('EXTRA',p.name)
  for l in gzip.open(p,'rt'):
   if l.startswith(('!Series_title','!Series_overall','!Sample_title','!Sample_characteristics')):print(l.strip())

