import pathlib,xml.etree.ElementTree as ET,zipfile,re
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2'
for fn in ['GC2025_stroma_myeloid.xml','GC2026_LysoPS.xml']:
 root=ET.parse(R/'01_metadata'/fn).getroot();txt=' '.join(root.itertext());print(fn,sorted(set(re.findall(r'(?:GSE|HRA|CRA|PXD|PRJCA)\d+',txt))))
 for e in root.iter('sec'):
  title=e.find('title')
  if title is not None and ('Data' in ''.join(title.itertext()) or 'data' in ''.join(title.itertext())):print(' '.join(e.itertext())[:2500])
 for e in root.iter('supplementary-material'):print('SUPP',ET.tostring(e,encoding='unicode')[:1200])
p=R/'01_metadata'/'GSE228598_supplement.zip'
with zipfile.ZipFile(p) as z:print('EXTERNAL_SUPP',[(x.filename,x.file_size) for x in z.infolist()])
