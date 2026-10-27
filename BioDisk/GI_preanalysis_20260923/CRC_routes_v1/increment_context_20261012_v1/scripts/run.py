from pathlib import Path
import subprocess,json,datetime,sys,traceback
R=Path(__file__).resolve().parent.parent
previous=R/'PROCESS_STATUS.json'
if previous.exists():previous.rename(R/('PROCESS_STATUS_'+datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S')+'.json'))
record={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':{}}
try:
 for name in ['06_report.py']:
  log=R/'logs'/('stage_'+name+'.log')
  if log.exists():log.rename(log.with_name(log.name+'.'+datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S')+'.previous'))
  with log.open('w') as f:
   p=subprocess.run([sys.executable,str(R/'scripts'/name)],stdout=f,stderr=subprocess.STDOUT)
  record['stages'][name]=p.returncode
  (R/'PROCESS_STATUS.json').write_text(json.dumps(record,indent=2))
  print(name,'EXIT',p.returncode,flush=True)
  if p.returncode:raise RuntimeError('stage failed '+name)
 record['all_success']=True
except Exception as e:record.update(all_success=False,error=repr(e));traceback.print_exc()
finally:
 record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 (R/'PROCESS_STATUS.json').write_text(json.dumps(record,indent=2))
sys.exit(0 if record['all_success'] else 1)
