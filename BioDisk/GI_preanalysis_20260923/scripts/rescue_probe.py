import pathlib,pandas as pd,gzip,re
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/01_metadata')
print((r/'GSE314812.files.html').read_text())
for l in gzip.open(r/'GSE314812.soft.gz','rt'):
 if l.startswith(('!Series_title','!Series_overall','!Sample_title','!Sample_characteristics','!Sample_supplementary','!Series_supplementary')):print(l.strip())
p=r/'GC2026_clinical.xlsx';e=pd.ExcelFile(p)
for sh in e.sheet_names:print('CLINICAL',sh,pd.read_excel(p,sheet_name=sh,header=None).head(10).to_string(index=False,header=False))
for fn in ['iprox_common.js','iprox_MSV02201.js','iprox_MSV02202.js']:
 t=(r/fn).read_text();print('JS',fn)
 for l in t.splitlines():
  if 'jsonServer =' in l or ('url:' in l and any(k in l.lower() for k in ['project','sample','protocol'])):print(l.strip())
