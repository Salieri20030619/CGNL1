from pathlib import Path
import subprocess,json,textwrap
R=Path(__file__).resolve().parent
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/SCP2079_qualification_20261014_v1'
receipt=json.loads((R/'BUNDLE_RECEIPT.json').read_text(encoding='utf-8'))
flags=getattr(subprocess,'CREATE_NO_WINDOW',0)
p=subprocess.run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20',str(R/receipt['path']),'sgs10:'+remote+'/'],capture_output=True,creationflags=flags)
if p.returncode:raise RuntimeError(p.stderr.decode('utf-8'))
code='''
from pathlib import Path
import json,zipfile,hashlib,datetime
R=Path(REMOTE)
zp=R/BUNDLE
assert hashlib.sha256(zp.read_bytes()).hexdigest()==EXPECTED
with zipfile.ZipFile(zp) as z:
 assert z.testzip() is None
 for name in z.namelist():
  dst=(R/name).resolve()
  assert R.resolve() in dst.parents
  if dst.exists():assert dst.read_bytes()==z.read(name)
 z.extractall(R)
manifest=json.loads((R/'FILE_MANIFEST.json').read_text(encoding='utf-8'))
for x in manifest:
 p=R/x['path']
 assert p.stat().st_size==x['bytes']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256']
status=json.loads((R/'PROGRESS_SNAPSHOT.json').read_text(encoding='utf-8'))
progress=R.parent/'CGNL1_PROGRESS_20261014_SCP2079.json'
payload=json.dumps(status,indent=2)
if progress.exists():assert progress.read_text(encoding='utf-8')==payload
else:progress.write_text(payload,encoding='utf-8')
receipt={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bundle_sha256':EXPECTED,'manifest_files_verified':len(manifest),'all_pass':True,'directory':str(R),'progress_file':str(progress),'matrix_pending_authenticated_download':True}
(R/'SERVER_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt))
'''
code='REMOTE='+repr(remote)+'\nBUNDLE='+repr(receipt['path'])+'\nEXPECTED='+repr(receipt['sha256'])+'\n'+code
p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3 -'],input=code.encode('utf-8'),capture_output=True,creationflags=flags)
if p.returncode:raise RuntimeError(p.stderr.decode('utf-8')+p.stdout.decode('utf-8'))
data=json.loads(p.stdout.decode('utf-8'))
(R/'SERVER_RECEIPT.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
print(json.dumps(data,indent=2),flush=True)