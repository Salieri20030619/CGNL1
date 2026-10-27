from pathlib import Path
import json,sys,subprocess,datetime,traceback
p=Path(sys.argv[1]);m=json.loads(p.read_text(encoding='utf-8'));record=dict(m)
try:
 with Path(m['log']).open('wb') as f:
  proc=subprocess.Popen([m['python'],m['script']],cwd=m['cwd'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,close_fds=True);record['worker_pid']=proc.pid;Path(m['state']).write_text(json.dumps(record,indent=2),encoding='utf-8');record['exit_code']=proc.wait()
except Exception:
 record['launcher_error']=traceback.format_exc();record['exit_code']=1
record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();Path(m['receipt']).write_text(json.dumps(record,indent=2),encoding='utf-8')
