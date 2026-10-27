
from pathlib import Path
import json,time,datetime,os
p=Path.cwd();meta=p/'metadata'
def read(f):
 try:return json.loads(Path(f).read_text(encoding='utf-8'))
 except (FileNotFoundError,json.JSONDecodeError):return None
while True:
 reg=read(meta/'PROJECT_JOB_REGISTRY_216.json')
 if not reg:time.sleep(900);continue
 if reg.get('stop_requested') or reg.get('project_goal_complete'):break
 jobs=[]
 for f in reg['jobs']:
  r=read(f)
  if not r:jobs.append({'record':f,'state':'missing_launch'});continue
  receipt=read(r['receipt']);live=Path('/proc',str(r['wrapper_pid'])).exists()
  jobs.append({'record':f,'state':'exited' if receipt is not None else ('running_or_waiting' if live else 'missing_exit_receipt'),'exit_code':receipt.get('exit_code') if receipt is not None else None,'receipt':r['receipt'],'log':r['log']})
 s={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'monitor_pid':os.getpid(),'interval_seconds':900,'project_goal_complete':False,'scientific_status':reg['scientific_status'],'jobs':jobs,'awaiting_scientific_decision':bool(jobs) and all(j['state']=='exited' for j in jobs),'needs_review':any(j['state']=='missing_exit_receipt' or (j['state']=='exited' and j['exit_code']!=0) for j in jobs)}
 tmp=meta/'PROJECT_PROGRESS_LATEST_216.json.writing';tmp.write_text(json.dumps(s,indent=2),encoding='utf-8');tmp.replace(meta/'PROJECT_PROGRESS_LATEST_216.json')
 with (meta/'PROJECT_PROGRESS_HISTORY_216.jsonl').open('a',encoding='utf-8') as h:h.write(json.dumps(s)+chr(10))
 print(json.dumps(s),flush=True);time.sleep(900)
