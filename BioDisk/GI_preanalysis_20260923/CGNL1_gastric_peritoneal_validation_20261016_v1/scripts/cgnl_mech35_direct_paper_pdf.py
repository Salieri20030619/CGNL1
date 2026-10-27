from pathlib import Path
import requests,subprocess,json,re,hashlib
O=Path.cwd()/'mechanism_investigation_v2';u='https://kumadai.repo.nii.ac.jp/record/2000526/files/igaku_kou2464zenbun.pdf';p=O/'literature/Yonemura2024_public_manuscript.pdf'
a=requests.get(u,timeout=(15,60));a.raise_for_status();assert a.content.startswith(b'%PDF');p.write_bytes(a.content);t=p.with_suffix('.txt');r=subprocess.run(['pdftotext','-layout',str(p),str(t)],capture_output=True);assert r.returncode==0;s=t.read_text(encoding='utf-8');ids=sorted(set(re.findall(r'GSE\d+|PRJ\w+\d+|DRA\d+|E-MTAB-\d+',s)));out={'url':u,'sha256':hashlib.sha256(a.content).hexdigest(),'bytes':len(a.content),'accessions':ids,'sections':[]}
for term in ['Data and code availability','RNA-seq data','Deposited data','CGNL1','mesothelial cells were treated']:
 for m in re.finditer(re.escape(term),s,re.I):out['sections'].append({'term':term,'text':s[max(0,m.start()-150):m.start()+1400]})
(O/'metadata/YONEMURA_PAPER_35.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
