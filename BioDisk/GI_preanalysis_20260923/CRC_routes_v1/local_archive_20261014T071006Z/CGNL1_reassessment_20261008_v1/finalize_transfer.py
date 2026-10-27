from pathlib import Path
import subprocess,time,json,hashlib,zipfile
root=Path(r'C:\Users\31937\CGNL1_reassessment_20261008_v1')
remote="""from pathlib import Path
import subprocess,json,hashlib,zipfile
r=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1')
with (r/'logs/report.log').open('w') as log:code=subprocess.call(['/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python',str(r/'build_reassessment_report.py')],stdout=log,stderr=subprocess.STDOUT)
(r/'logs/report.exit').write_text(str(code)+'\\n');assert code==0,(r/'logs/report.log').read_text()[-3000:]
files=[p for p in sorted(r.rglob('*')) if p.is_file() and p.name not in ['SHA256SUMS.txt','CGNL1_reassessment_bundle.zip']]
(r/'SHA256SUMS.txt').write_text('\\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(r).as_posix() for p in files)+'\\n')
with zipfile.ZipFile(r/'CGNL1_reassessment_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in files+[r/'SHA256SUMS.txt']:z.write(p,p.relative_to(r))
with zipfile.ZipFile(r/'CGNL1_reassessment_bundle.zip') as z:assert z.testzip() is None
p=r/'CGNL1_reassessment_bundle.zip';print(json.dumps({'files':len(files),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}))
"""
def run_retry(args,payload=None):
 for attempt in range(1,5):
  result=subprocess.run(args,input=payload,stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=subprocess.CREATE_NO_WINDOW)
  if result.returncode==0:return result.stdout.decode('utf-8',errors='replace')
  print(json.dumps({'attempt':attempt,'exit':result.returncode,'stderr':result.stderr.decode('utf-8',errors='replace')[-700:]}),flush=True)
  if attempt<4:time.sleep(15)
 raise RuntimeError('Remote transfer did not complete after four controlled attempts')
output=run_retry(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','sgs10','python3','-'],remote.encode('utf-8'));expected=json.loads(output.strip().splitlines()[-1]);print(json.dumps(expected),flush=True)
zpath=root/'CGNL1_reassessment_bundle.zip';run_retry(['scp','-q','sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1/CGNL1_reassessment_bundle.zip',str(zpath)])
assert hashlib.sha256(zpath.read_bytes()).hexdigest()==expected['sha256']
dest=(root/'delivery').resolve()
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for name in z.namelist():assert (dest/name).resolve().is_relative_to(dest)
 z.extractall(dest)
count=0
for line in (dest/'SHA256SUMS.txt').read_text().splitlines():
 digest,rel=line.split('  ',1);assert hashlib.sha256((dest/rel).read_bytes()).hexdigest()==digest,rel;count+=1
record={'files_verified':count,'failed':0,'zip_sha256':expected['sha256']};(root/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(record,indent=2));print(json.dumps(record),flush=True)
