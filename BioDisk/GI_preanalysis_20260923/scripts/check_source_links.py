import pathlib,xml.etree.ElementTree as ET
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
t=ET.parse(r/'01_metadata'/'GC_ascites_2023.xml').getroot()
for e in t.iter():
 if e.tag in ['ext-link','supplementary-material','media'] and ('Source' in ''.join(e.itertext()) or 'source' in ''.join(e.itertext()) or 'MOESM' in str(e.attrib)):
  print(ET.tostring(e,encoding='unicode')[:2500])
