from pathlib import Path
import json,hashlib,datetime,zipfile,base64,io
from bs4 import BeautifulSoup
from PIL import Image
O=Path.cwd()
def save(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
f=O/'report/CGNL1_gastric_peritoneal_validation.html';b=f.read_bytes();s=b.decode('utf-8');assert not b.startswith(b'\xef\xbb\xbf') and '\ufffd' not in s
soup=BeautifulSoup(s,'html.parser');imgs=soup.find_all('img');assert len(imgs)==2
for im in imgs: Image.open(io.BytesIO(base64.b64decode(im['src'].split(',',1)[1]))).verify()
assert json.loads((O/'NUMERICAL_QA.json').read_text())['all_pass']
save(O/'REPORT_QA.json',{'utf8_no_bom':True,'replacement_characters':0,'embedded_images_verified':2,'tables':len(soup.find_all('table')),'visual_review':'Chinese HTML opening view and both complete figures inspected; readable','html_sha256':sha(f)})
status={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'actual_data_validation_complete_pending_user_review','whole_project_complete':False,'server_root':str(O),'report':str(f),'decision':json.loads((O/'DECISION.json').read_text(encoding='utf-8'))}
save(O/'STATUS.json',status)
save(O.parent/'CGNL1_CURRENT_STATE_20261016_GASTRIC_PERITONEAL.json',status)
files=[]
for folder in ['report','figures','results','scripts','logs']:
 files.extend(p for p in (O/folder).rglob('*') if p.is_file())
files.extend(O.glob('*.json'))
files.extend(p for p in (O/'metadata').glob('*') if p.is_file() and p.suffix in ['.json','.tsv','.txt'])
files.extend([O/'processed/original31_sample_values.tsv',O/'processed/bulk_sample_measurements.tsv'])
files=sorted(set(files))
manifest={str(p.relative_to(O)): {'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}
save(O/'DELIVERY_MANIFEST.json',manifest)
zpath=O/'report/CGNL1_gastric_peritoneal_evidence_bundle.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
 z.write(O/'DELIVERY_MANIFEST.json','DELIVERY_MANIFEST.json')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'zip_sha256':sha(zpath),'html_sha256':sha(f),'zip_bytes':zpath.stat().st_size,'manifest_entries':len(manifest),'archive_crc_pass':True,'raw_and_intermediate_data_location':str(O),'scope':'This round complete; final scientific mainline not established; user review pending'}
save(O/'VALIDATION_COMPLETE.json',receipt);print(json.dumps(receipt,ensure_ascii=False))
