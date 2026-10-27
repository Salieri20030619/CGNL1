from pathlib import Path
import zipfile,subprocess,shutil,json,datetime,re
O=Path.cwd();z=O.parent/'GC_validation_v2/01_metadata/GSE228598_supplement.zip';p=O/'metadata/GSE228598_supplement.pdf'
with zipfile.ZipFile(z) as f:entries=[i for i in f.infolist() if i.filename.endswith('.pdf')];assert len(entries)==1;p.write_bytes(f.read(entries[0]))
cmd=shutil.which('pdftotext');render=shutil.which('pdftoppm');info=shutil.which('pdfinfo');record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_archive':str(z),'pdf':str(p),'pdftotext':cmd,'pdftoppm':render}
if cmd:
 r=subprocess.run([cmd,'-layout',str(p),'-'],capture_output=True,timeout=30);assert r.returncode==0;txt=r.stdout.decode('utf-8');(O/'metadata/GSE228598_supplement_text.txt').write_text(txt,encoding='utf-8');record['text_pages']=len(txt.split('\f'));print(txt[:15000])
if render:
 r=subprocess.run([render,'-f','1','-singlefile','-scale-to','1500','-png',str(p),str(O/'metadata/GSE228598_supplement_page1')],capture_output=True,timeout=30);record['render_exit']=r.returncode
(O/'metadata/PDF_INSPECTION.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record))
