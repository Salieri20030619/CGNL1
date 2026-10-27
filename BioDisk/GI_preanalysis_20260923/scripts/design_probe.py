import pathlib,gzip,json
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for p in list((r/'01_metadata').glob('GSE*.soft.gz'))+[r/'02_raw'/'GSE200697_family.soft.gz']:
 if not p.exists():continue
 print('\nDATA',p.name)
 with gzip.open(p,'rt') as f:
  table=False;c=0
  for l in f:
   if l.startswith(('!Series_title','!Series_overall','^SAMPLE','!Sample_title','!Sample_characteristics')):print(l.strip())
   if l.startswith('!platform_table_begin') or l.startswith('!sample_table_begin'):table=True;c=0
   elif table and c<2:print('TABLE',l.strip()[:1500]);c+=1
   elif l.startswith('!platform_table_end') or l.startswith('!sample_table_end'):table=False
