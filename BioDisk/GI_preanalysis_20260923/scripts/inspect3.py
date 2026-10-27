import pathlib,gzip,re,xml.etree.ElementTree as ET,json
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/01_metadata')
for p in r.glob('*.soft.gz'):
 print('\nDATASET',p.name)
 for line in gzip.open(p,'rt'):
  if line.startswith(('!Series_title','!Series_overall','!Series_pubmed','^SAMPLE','!Sample_title','!Sample_characteristics')):print(line.strip())
for p in r.glob('*.files.html'):
 if p.stem.startswith(('GSE163','GSE234','GSE200','GSE228')):print(p.name,p.read_text())
p=r/'PXD064711_central.xml'
if p.exists():print(p.read_text()[:10000])
p=r/'GSE145356_counts.txt.gz'
if p.exists():
 with gzip.open(p,'rt') as f: print('COUNT HEAD',[next(f) for _ in range(3)])
p=r/'GSE294559.soft.txt'
if p.exists():
 for line in p.read_text().splitlines():
  if line.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics')):print(line)
