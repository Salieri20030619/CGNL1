from pathlib import Path
import subprocess,zipfile,json,hashlib
L=Path(__file__).resolve().parent
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/increment_context_20261012_v1'
for name in ['CGNL1_increment_context_bundle.zip','BUNDLE_RECEIPT.json']:
 subprocess.run(['scp','-q','sgs10:'+remote+'/'+name,str(L/name)],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
receipt=json.loads((L/'BUNDLE_RECEIPT.json').read_text())
bundle=L/'CGNL1_increment_context_bundle.zip'
assert hashlib.sha256(bundle.read_bytes()).hexdigest()==receipt['sha256']
D=L/'delivery';D.mkdir(exist_ok=True)
with zipfile.ZipFile(bundle) as z:
 for member in z.infolist():
  p=(D/member.filename).resolve();assert p.is_relative_to(D.resolve())
 z.extractall(D)
checked=0
for line in (D/'SHA256SUMS.txt').read_text().splitlines():
 h,p=line.split('  ',1);assert hashlib.sha256((D/p).read_bytes()).hexdigest()==h,p;checked+=1
check=dict(bundle_sha256_match=True,manifest_files_verified=checked,report=str(D/'CGNL1_increment_context_report.html'))
(L/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(check,indent=2))
print(json.dumps(check))
