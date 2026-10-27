from pathlib import Path
import subprocess,sys,json,datetime,os,traceback
R=Path(__file__).resolve().parent.parent
(R/'logs').mkdir(exist_ok=True)
state={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'host':'sgs10/osv01','python':sys.executable,'stages':[]}
(R/'JOB_STATUS.json').write_text(json.dumps(state,indent=2),encoding='utf-8')
rc=0
try:
 for stage in ['01_matrix_qualify.py','02_freeze_extract.py','03_analyze.py','04_qa.py']:
  info={'script':stage,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  with (R/'logs'/(stage+'.log')).open('w',encoding='utf-8') as log:
   p=subprocess.run([sys.executable,'-u',str(R/'scripts'/stage)],stdout=log,stderr=subprocess.STDOUT)
  info.update(exit_code=p.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  (R/'logs'/(stage+'.exit.json')).write_text(json.dumps(info,indent=2),encoding='utf-8');state['stages'].append(info)
  (R/'JOB_STATUS.json').write_text(json.dumps(state,indent=2),encoding='utf-8')
  if p.returncode:rc=p.returncode;break
except Exception:
 state['error']=traceback.format_exc();rc=1
state.update(exit_code=rc,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'JOB_COMPLETE.json').write_text(json.dumps(state,indent=2),encoding='utf-8')
sys.exit(rc)