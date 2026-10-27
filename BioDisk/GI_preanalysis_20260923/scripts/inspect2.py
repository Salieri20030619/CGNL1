import pathlib,xml.etree.ElementTree as ET,importlib.util,json
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for f in ['GC_ascites_2023.xml','CRC_multiome_2025.xml']:
 root=ET.parse(r/'01_metadata'/f).getroot()
 for e in root.iter():
  if e.tag in ['p','sec'] and ('HRA002712' in ''.join(e.itertext()) or (e.tag=='p' and 'nine' in ''.join(e.itertext()))): print(f,' '.join(e.itertext())[:3500])
 for e in root.iter('supplementary-material'):print(f,ET.tostring(e,encoding='unicode')[:1000])
p=r/'01_metadata'/'GSE145356.soft.txt'
if p.exists():
 for line in p.read_text().splitlines():
  if line.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics','!Series_pubmed','!Series_overall')):print(line)
print('PACKAGES',[(x,importlib.util.find_spec(x) is not None) for x in ['numpy','pandas','scipy','h5py','anndata','scanpy','matplotlib']])
