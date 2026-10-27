from pathlib import Path
import requests,zipfile,subprocess,json,re
O=Path.cwd()/'mechanism_investigation_v2';u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12013559/supplementaryFiles';a=requests.get(u,timeout=(15,90));a.raise_for_status();p=O/'metadata/GSE251950_paper_supplement.zip';p.write_bytes(a.content);out=[]
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None;print('FILES',[(i.filename,i.file_size) for i in z.infolist()])
 for i in z.infolist():
  if i.filename.endswith('.pdf'):
   f=O/'literature'/Path(i.filename).name;f.write_bytes(z.read(i.filename));t=f.with_suffix('.txt');subprocess.run(['pdftotext','-layout',str(f),str(t)],check=True);s=t.read_text(encoding='utf-8');hits=[]
   for n,page in enumerate(s.split('\f'),1):
    for m in re.finditer(r'GC6|GC 6|peritone|metastatic site|metastasis site',page,re.I):hits.append({'page':n,'context':page[max(0,m.start()-180):m.start()+500]})
   out.append({'pdf':str(f),'hits':hits});print('HITS',json.dumps(out[-1],ensure_ascii=False)[:4500])
(O/'metadata/SPATIAL_SUPPLEMENT_87.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
