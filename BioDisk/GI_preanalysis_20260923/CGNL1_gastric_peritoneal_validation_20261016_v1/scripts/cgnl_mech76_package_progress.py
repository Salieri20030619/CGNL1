from pathlib import Path
import json,hashlib,zipfile,datetime
O=Path.cwd()/'mechanism_investigation_v2';files=[O/'report/CGNL1_mechanism_progress_20261016.html',O/'report/ascites_identity_comparison.png'];files+=list((O/'results').glob('*'));files+=list(O.glob('*SPECIFICATION*.json'));files +=[O/'FPKM_RECONSTRUCTION_AMENDMENT.json',O/'STATUS.json'];files=[p for p in sorted(set(files)) if p.is_file()];manifest={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
p=O/'report/CGNL1_mechanism_progress_evidence.zip'
with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,str(f.relative_to(O)))
 z.writestr('MANIFEST.json',json.dumps(manifest,indent=2).encode('utf-8'))
with zipfile.ZipFile(p) as z:assert z.testzip() is None
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'zip_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'html_sha256':manifest['report/CGNL1_mechanism_progress_20261016.html'],'files':len(files),'bytes':p.stat().st_size,'status':'interim evidence update; mechanism goal remains active'};(O/'report/PROGRESS_DELIVERY_RECEIPT.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
