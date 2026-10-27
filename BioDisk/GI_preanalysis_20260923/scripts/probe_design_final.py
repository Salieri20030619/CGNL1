import pathlib,gzip,pandas as pd
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for fn in ['GSE163171.soft.gz','GSE187650.soft.gz','GSE294559.soft.gz']:
 p=r/'01_metadata'/fn
 if p.exists():
  print('DATA',fn)
  for l in gzip.open(p,'rt'):
   if l.startswith(('!Series_title','!Series_overall','!Sample_title')):print(l.strip())
for fn in ['GSE163171.files.html','GSE187650.files.html']:
 p=r/'01_metadata'/fn
 if p.exists():print(p.read_text())
p=r/'02_raw'/'GSE200697_family.soft.gz'
with gzip.open(p,'rt') as f:
 for l in f:
  if l.startswith(('!Sample_data_processing','!Sample_supplementary_file')):print(l.strip())
