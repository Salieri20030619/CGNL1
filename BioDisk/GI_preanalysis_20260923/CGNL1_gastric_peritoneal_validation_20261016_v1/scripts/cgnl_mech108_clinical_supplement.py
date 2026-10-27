from pathlib import Path
import json,zipfile,io,hashlib,subprocess
import openpyxl
O=Path.cwd()/'mechanism_investigation_v2'
p=O/'metadata/GSE314812_article_supplement.zip'
z=zipfile.ZipFile(p);assert z.testzip() is None
out={'archive_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'members':z.namelist(),'matches':[],'sheet_overviews':[]}
for name in z.namelist():
 if not name.endswith('.xlsx'):continue
 wb=openpyxl.load_workbook(io.BytesIO(z.read(name)),read_only=True,data_only=True)
 for ws in wb:
  rows=list(ws.iter_rows(values_only=True))
  out['sheet_overviews'].append({'file':name,'sheet':ws.title,'rows':len(rows),'first_rows':rows[:4]})
  for i,row in enumerate(rows):
   if any('GC05' in str(x) or 'GC5' in str(x) for x in row):out['matches'].append({'file':name,'sheet':ws.title,'row':i+1,'context':rows[max(0,i-2):i+3]})
 wb.close()
(O/'metadata/GSE314812_CLINICAL_SUPPLEMENT_108.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
print('CLINICAL',json.dumps(out,ensure_ascii=False,default=str)[:24000])
qual=json.loads((O/'metadata/GSE308231_SPECIES_AND_METADATA_106.json').read_text());print('SPECIES',json.dumps(qual['feature_qualification']))
for name in ['INDEX_COMPLETE.json','DRA_QUANTIFICATION_COMPLETE.json','GSE308231_DOWNLOAD.json','SOURCE_COMPLETE.json','DRA_DOWNLOAD_COMPLETE.json']:
 for p in O.rglob(name):print('MARKER',str(p),p.read_text()[:4000])
print('JOBS',subprocess.run(['ps','-u','wangmingyuan','-o','pid,etime,pcpu,pmem,args'],capture_output=True,text=True).stdout[-25000:])
