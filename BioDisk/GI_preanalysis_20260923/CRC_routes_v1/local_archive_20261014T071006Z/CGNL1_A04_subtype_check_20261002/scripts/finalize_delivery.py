
from pathlib import Path
import json,hashlib,tarfile,zipfile,subprocess,time,shutil
r=Path(r'C:\Users\31937\CGNL1_A04_subtype_check_20261002')
files=[r/'A04检查说明.txt',r/'MODULE_STATUS.json',r/'scripts/build_report.py',r/'figures/A04_checks.png',r/'figures/A04_checks.svg',r/'results/fixed_baseline_subtype_diagnostic.tsv']
(r/'DELIVERY_QA.json').write_text(json.dumps({'analysis_exit_codes':[0,0],'previous_profile_reproduction':'passed','identity_freeze_sha256':'passed','download_44_hashes':'passed','figures_visually_reviewed':True,'all_zero_unestimable_model_checked':True,'scope':'two-cohort diagnostic, not full project completion'},indent=2));files.append(r/'DELIVERY_QA.json')
a={p.relative_to(r).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files};(r/'ADDENDUM_MANIFEST.json').write_text(json.dumps(a,indent=2));files.append(r/'ADDENDUM_MANIFEST.json')
with tarfile.open(r/'report_addendum.tar.gz','w:gz') as tf:
 for p in files:tf.add(p,arcname=p.relative_to(r).as_posix())
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/A04_subtype_check_20261002_v1'
for attempt in range(3):
 with (r/'report_addendum.tar.gz').open('rb') as inp:x=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','cat > '+remote+'/report_addendum.tar.gz'],stdin=inp,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode==0:break
 time.sleep(3)
assert x.returncode==0,x.stderr
code="from pathlib import Path\nimport tarfile,json,hashlib\nr=Path("+repr(remote)+")\nwith tarfile.open(r/'report_addendum.tar.gz') as tf:\n for m in tf.getmembers():assert (r/m.name).resolve().is_relative_to(r.resolve()) and not m.issym() and not m.islnk()\n tf.extractall(r)\na=json.loads((r/'ADDENDUM_MANIFEST.json').read_text());assert all(hashlib.sha256((r/f).read_bytes()).hexdigest()==h for f,h in a.items());print('REMOTE_VERIFIED',len(a))"
for attempt in range(3):
 x=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3 -'],input=code.encode(),capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode==0:break
 time.sleep(3)
assert x.returncode==0,x.stderr;print(x.stdout.decode())
files=[p for p in r.rglob('*') if p.is_file() and not p.name.endswith('.tar.gz') and p.name!='FINAL_MANIFEST.json'];a={p.relative_to(r).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files};(r/'FINAL_MANIFEST.json').write_text(json.dumps(a,indent=2));files.append(r/'FINAL_MANIFEST.json')
out=Path(r'C:\Users\31937\Downloads\CGNL1_A04_subtype_check_20261002.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,arcname=p.relative_to(r).as_posix())
with zipfile.ZipFile(out) as z:assert z.testzip() is None
shutil.copy2(r/'A04检查说明.txt',Path(r'C:\Users\31937\Downloads\CGNL1_A04检查说明_20261002.txt'))
print('DELIVERED',out,'FILES',len(files))
