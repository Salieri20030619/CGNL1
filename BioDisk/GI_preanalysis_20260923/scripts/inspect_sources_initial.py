import pathlib,json,re,xml.etree.ElementTree as ET,importlib.util
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for p in (r/'01_metadata').glob('*'):
 if p.suffix=='.xml':
  root=ET.parse(p).getroot(); text=' '.join(root.itertext()); print(p.name,'ACCESSIONS',sorted(set(re.findall(r'(?:GSE|HRA|CRA|PRJCA|EGAS)\d+',text))))
  for e in root.iter('sec'):
   ti=e.find('title')
   if ti is not None and 'data availab' in ''.join(ti.itertext()).lower():print(' '.join(e.itertext()))
 elif p.name.endswith('search.json'): print(p.name,p.read_text())
 elif 'zenodo' in p.name:
  j=json.loads(p.read_text());print('ZENODO',[(f['key'],f['size'],f['links']) for f in j['files']])
 elif 'files.html' in p.name:print(p.name,p.read_text())
print('PACKAGES',[(x,importlib.util.find_spec(x) is not None) for x in ['numpy','pandas','scipy','h5py','anndata','scanpy','matplotlib']])
