from pathlib import Path
import subprocess,json,time,sys
L=Path(__file__).resolve().parent
R='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1'
files={name:(L/name).read_text(encoding='utf-8') for name in ['transfer_diagnostic.py','07_report.py','08_package.py']}
bootstrap="from pathlib import Path\nimport json,subprocess,sys\nR=Path("+repr(R)+")\nfiles=json.loads("+repr(json.dumps(files))+")\nfor name,code in files.items():(R/'scripts'/name).write_text(code,encoding='utf-8')\nPY='/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python'\nfor name in ['07_report.py']:\n p=subprocess.run([PY,str(R/'scripts/run.py'),name]);assert p.returncode==0,name\np=subprocess.run([PY,str(R/'scripts/08_package.py')]);assert p.returncode==0\nprint((R/'BUNDLE_RECEIPT.json').read_text())\n"
for attempt in range(3):
 p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3','-'],input=bootstrap.encode(),capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 sys.stdout.buffer.write(p.stdout);sys.stderr.buffer.write(p.stderr);sys.stdout.flush()
 if p.returncode==0:break
 if p.returncode!=255:p.check_returncode()
 if attempt<2:time.sleep(30)
p.check_returncode()
