from pathlib import Path
import json,subprocess,hashlib,datetime,re
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();rows=json.loads((O/'metadata/ORIGINAL31_FILE_LIST.json').read_text());files=[]
for r in rows:
 for f in r['files']:files.append(dict(f,folder_name=r['folder']['name']))
(O/'raw/original31').mkdir(exist_ok=True)
def download(f):
 p=O/'raw/original31'/f['filename'];u=f['content_details']['download_url'];r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','100','--retry','1','-sS','-w','%{http_code}','-o',str(p),u],capture_output=True,timeout=210);b=p.read_bytes() if p.exists() else b'';h=hashlib.sha256(b).hexdigest();expected=f['content_details']['sha256_hash'];rec={'filename':f['filename'],'folder_name':f['folder_name'],'url':u,'exit':r.returncode,'http':r.stdout.decode('utf-8'),'bytes':len(b),'sha256':h,'expected_sha256':expected,'hash_matches':h==expected,'magic_hex':b[:12].hex()};print('CHP',rec['filename'],rec['exit'],rec['hash_matches'],flush=True);return rec
with ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(download,files))
(O/'metadata/ORIGINAL31_DOWNLOAD_RECEIPTS.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8');assert all(x['hash_matches'] for x in records)
note={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':len(records),'tumor_with_PC':sum(x['filename'].startswith('Tumor with PC ') for x in records),'tumor_without_PC':sum(x['filename'].startswith('Tumor without PC ') for x in records),'healthy':sum(x['filename'].startswith('Healthy') for x in records),'label_conflict':'folder named Tumors without peritoneal carcinomatosis contains12filesnamedTumor with PC;otherfolderNew foldercontains19filesnamedTumor without PC','resolution_pending':'inspect file-internal metadata before assigning clinical groups; do not select label interpretation using target expression','healthy_pairing':'no verified patient-ID matching; do not assume23healthy samplesarepairedto31tumors'}
(O/'metadata/ORIGINAL31_QUALIFICATION.json').write_text(json.dumps(note,indent=2),encoding='utf-8');print(json.dumps(note));
sample=next((O/'raw/original31').glob('Tumor with PC *.chp'));b=sample.read_bytes();strings=re.findall(rb'[\x20-\x7e]{6,}',b[:15000]);print('ASCII_HEADER',[s.decode('ascii') for s in strings[:45]]);print('MAGIC',b[:30].hex())
