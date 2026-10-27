import pathlib,zipfile,json,pypdf
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');M=O/'metadata'
z=zipfile.ZipFile(M/'GSE225857_supplements.zip');print('ZIP',z.namelist())
for n in z.namelist():
 if n.endswith('.pdf'):
  p=M/pathlib.Path(n).name;p.write_bytes(z.read(n));r=pypdf.PdfReader(p)
  for i,page in enumerate(r.pages):
   t=page.extract_text() or ''
   if 'Table S1' in t or 's1231' in t or 'Patient' in t:print('PAGE',i+1,t[:6000])
j=json.loads((M/'spatial_table2.json').read_text(encoding='utf-8'));print('CLINICAL_SPATIAL_TAB1');print(j['Tab 1'])
