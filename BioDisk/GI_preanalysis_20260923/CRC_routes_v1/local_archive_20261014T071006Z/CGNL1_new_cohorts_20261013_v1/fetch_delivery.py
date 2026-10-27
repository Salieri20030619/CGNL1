from pathlib import Path
import subprocess,zipfile,json,hashlib
L=Path(__file__).resolve().parent
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1'
def run(a):
 for attempt in range(3):
  p=subprocess.run(a,capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if p.returncode==0:return
  print(p.stderr[-1000:],flush=True)
 p.check_returncode()
for name in ['CGNL1_new_cohorts_bundle.zip','BUNDLE_RECEIPT.json']:
 run(['scp','-q','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10:'+remote+'/'+name,str(L/name)])
rec=json.loads((L/'BUNDLE_RECEIPT.json').read_text());b=L/'CGNL1_new_cohorts_bundle.zip'
assert hashlib.sha256(b.read_bytes()).hexdigest()==rec['sha256']
D=L/'delivery';D.mkdir(exist_ok=True)
with zipfile.ZipFile(b) as z:
 for m in z.infolist():assert (D/m.filename).resolve().is_relative_to(D.resolve())
 z.extractall(D)
n=0
for line in (D/'SHA256SUMS.txt').read_text().splitlines():
 h,p=line.split('  ',1);assert hashlib.sha256((D/p).read_bytes()).hexdigest()==h,p;n+=1
out=dict(bundle_sha256_match=True,manifest_files_verified=n,report=str(D/'CGNL1_new_cohorts_report.html'))
(L/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
