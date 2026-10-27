from pathlib import Path
import subprocess,json
H=Path(r'D:\CGNL1\handoff')
ready=json.loads((H/'RAW_SERVER_READY.json').read_text(encoding='utf-8'));assert ready['all_pass']
code='''
from pathlib import Path
import subprocess,json,os
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/SCP2079_paired_20261014_v1')
assert not (R/'JOB_COMPLETE.json').exists() and not (R/'JOB_STATUS.json').exists()
for p in (R/'raw/expression_epi_counts').iterdir():p.chmod(0o444)
with (R/'logs/launcher.log').open('wb') as log:
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python','-u',str(R/'scripts/run_pipeline.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(R/'JOB_PID.txt').write_text(str(p.pid),encoding='utf-8')
print(json.dumps({'pid':p.pid,'directory':str(R),'host':'sgs10/osv01'}))
'''
p=subprocess.run(['ssh','-o','BatchMode=yes','sgs10','python3 -'],input=code.encode('utf-8'),capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
assert p.returncode==0,p.stderr.decode('utf-8')
r=json.loads(p.stdout.decode('utf-8'));(H/'ANALYSIS_LAUNCH.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True)