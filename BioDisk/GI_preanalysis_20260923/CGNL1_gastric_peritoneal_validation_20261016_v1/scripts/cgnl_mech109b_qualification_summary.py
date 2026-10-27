from pathlib import Path
import json,zipfile,io,re
import subprocess,openpyxl
O=Path.cwd()/'mechanism_investigation_v2';z=zipfile.ZipFile(O/'metadata/GSE314812_article_supplement.zip');hits=[]
for name in z.namelist():
 if not name.endswith('.pdf'):continue
 p=O/'literature'/name;p.write_bytes(z.read(name));s=subprocess.run(['pdftotext','-layout',str(p),'-'],capture_output=True,check=True).stdout.decode('utf-8')
 for i,text in enumerate(s.split('\f')):

  if any(x in text.lower() for x in ['gc05','longitudinal','repeated','18 pairs','twice']):hits.append({'file':name,'page':i+1,'text':text})
out={'pdf_hits':hits,'clinical_rows':[]}
for n in ['42003_2026_9906_MOESM3_ESM.xlsx','42003_2026_9906_MOESM17_ESM.xlsx']:
 wb=openpyxl.load_workbook(io.BytesIO(z.read(n)),read_only=True,data_only=True)
 for ws in wb:
  rows=list(ws.values);chosen=[{'row':i+1,'values':r} for i,r in enumerate(rows) if i<4 or 'GC05' in r]
  out['clinical_rows'].append({'file':n,'sheet':ws.title,'rows':chosen})
(O/'metadata/CLINICAL_QUALIFICATION_109.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str),encoding='utf-8');print(json.dumps(out,ensure_ascii=False,default=str)[:18000])
q=json.loads((O/'metadata/GSE308231_SPECIES_AND_METADATA_106.json').read_text());print('FEATURE',q['feature_qualification'])
for directory,pattern in [(O.parent/'logs','cgnl_mech67*.log'),(O.parent/'logs','cgnl_mech85*.log'),(O/'logs','salmon_index.log')]:
 for p in sorted(directory.glob(pattern))[-1:]:print('LOG',p.name,p.read_text()[-1800:])
