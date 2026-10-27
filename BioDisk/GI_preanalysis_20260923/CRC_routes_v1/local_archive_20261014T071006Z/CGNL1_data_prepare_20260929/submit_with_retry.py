from pathlib import Path
import base64,subprocess,time,json
r=Path(r'C:/Users/31937/CGNL1_data_prepare_20260929')
code=base64.b64encode((r/'download_manager.py').read_bytes()).decode()
remote="""import base64,subprocess,json,os
from pathlib import Path
r=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/data_prepare_20260929')
r.mkdir(exist_ok=True)
pidfile=r/'wave01.pid'
if pidfile.exists():
 try:
  os.kill(int(pidfile.read_text()),0);print('ALREADY_RUNNING',pidfile.read_text());raise SystemExit(0)
 except ProcessLookupError:pass
(r/'download_manager.py').write_bytes(base64.b64decode('"""+code+"""'))
launcher='import subprocess,json; from pathlib import Path; r=Path(__file__).parent; p=subprocess.run([\"python3\",str(r/\"download_manager.py\")]); (r/\"wave01.EXIT.json\").write_text(json.dumps({\"exit_code\":p.returncode}))'
(r/'launch_wave01.py').write_text(launcher)
with (r/'wave01.stdout.log').open('ab') as out,(r/'wave01.stderr.log').open('ab') as err:
 p=subprocess.Popen(['python3',str(r/'launch_wave01.py')],stdout=out,stderr=err,start_new_session=True,cwd=r)
pidfile.write_text(str(p.pid));print('SUBMITTED',p.pid)
"""
for attempt in range(1,11):
 p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','sgs10','python3','-'],input=remote,text=True,encoding='utf-8',capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 rec={'attempt':attempt,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'time':time.time()}
 (r/'submission_status.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(rec,ensure_ascii=False),flush=True)
 if p.returncode==0:break
 time.sleep(90)
else:raise SystemExit(1)
