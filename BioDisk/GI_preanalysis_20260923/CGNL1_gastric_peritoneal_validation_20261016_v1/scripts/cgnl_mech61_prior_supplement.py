from pathlib import Path
import zipfile,subprocess,json,re
O=Path.cwd()/'mechanism_investigation_v2';p=O/'metadata/DDR2_2016_supplement.zip'
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None;files=[{'name':i.filename,'bytes':i.file_size} for i in z.infolist()];print('FILES',json.dumps(files))
 for i in z.infolist():
  if i.filename.endswith('.pdf') or ('Fig1' in i.filename and i.filename.endswith('.jpg')):(O/'literature'/Path(i.filename).name).write_bytes(z.read(i.filename))
out=[]
for p in (O/'literature').glob('srep22371*.pdf'):
 t=p.with_suffix('.txt');subprocess.run(['pdftotext','-layout',str(p),str(t)],check=True);s=t.read_text(encoding='utf-8');pages=s.split('\f');hits=[]
 for n,page in enumerate(pages,1):
  for m in re.finditer('CGNL1',page,re.I):hits.append({'page':n,'context':page[max(0,m.start()-350):m.start()+350]})
 out.append({'pdf':str(p),'page_count':len(pages)-1,'target_hits':hits});print(json.dumps(out[-1],ensure_ascii=False))
(O/'metadata/DDR2_SUPPLEMENT_CHECK_61.json').write_text(json.dumps({'files':files,'PDFs':out},indent=2),encoding='utf-8')
