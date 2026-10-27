from pathlib import Path
import zipfile,json,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';p=O/'metadata/ascites_supplement_package.zip'
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 rows=[{'name':i.filename,'bytes':i.file_size} for i in z.infolist()];print(json.dumps(rows,indent=2))
 for i in z.infolist():
  if i.filename.endswith('.xlsx'):
   n=Path(i.filename).name;(O/'metadata'/n).write_bytes(z.read(i.filename))
(O/'metadata/ASCITES_SUPPLEMENT_INVENTORY.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
# Read-only OOXML inventory, no author spreadsheet modification.
for p in (O/'metadata').glob('41467_2023_36310*.xlsx'):
 with zipfile.ZipFile(p) as z:
  assert z.testzip() is None
  root=ET.fromstring(z.read('xl/workbook.xml'));print('WORKBOOK',p.name,[(e.attrib) for e in root.iter() if e.tag.endswith('}sheet')])
  for n in z.namelist():
   if n.startswith('xl/worksheets/sheet') and n.endswith('.xml'):
    root=ET.fromstring(z.read(n));dims=[e.attrib for e in root.iter() if e.tag.endswith('}dimension')];print('SHEET',n,dims)
