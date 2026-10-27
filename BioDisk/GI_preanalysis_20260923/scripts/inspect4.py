import pathlib,gzip,json,xml.etree.ElementTree as ET,pandas as pd
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for p in (r/'01_metadata').glob('*.xml'):
 if p.name in ['CRC_diapause.xml','GC_2021.xml']:
  rt=ET.parse(p).getroot()
  for e in rt.iter():
   if e.tag=='p' and any(x in ''.join(e.itertext()) for x in ['GSE','Data and Code Availability','GSA']):print(p.name,' '.join(e.itertext())[:1800])
for p in (r/'02_raw').glob('*'):
 print('FILE',p.name,p.stat().st_size)
 if p.suffix=='.gz' and any(k in p.name for k in ['meta','counts']):
  with gzip.open(p,'rt') as f:print('HEAD',[next(f).strip() for _ in range(3)])
p=r/'01_metadata'/'GC_ascites_source.xlsx'
if p.exists():
 x=pd.ExcelFile(p);print('SOURCE_SHEETS',x.sheet_names)
for p in (r/'logs').glob('*acquisition*'):print(p.name)
