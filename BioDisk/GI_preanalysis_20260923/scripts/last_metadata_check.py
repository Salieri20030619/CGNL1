import pathlib,json,gzip,re,html,xml.etree.ElementTree as ET
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
p=r/'01_metadata'/'iprox_proxi.json'
if p.exists():print('IPROX',p.read_text()[:14000])
p=r/'01_metadata'/'iprox_project.html'
if p.exists():
 s=p.read_text();print('IPROX_SCRIPT_URLS',re.findall(r'<script[^>]*src=[\"\x27]([^\"\x27]+)',s));print('IPROX_DESIGN',[(k,s.find(k)) for k in ['sample','Sample','protocol','Treatment','AGPAT2']])
for fn in ['GSE163171_family.soft.gz','GSE187650_family.soft.gz']:
 p=r/'01_metadata'/fn
 if p.exists():
  print('EXTRA',fn)
  for l in gzip.open(p,'rt'):
   if l.startswith(('!Series_title','!Series_overall','!Sample_title')):print(l.strip())
