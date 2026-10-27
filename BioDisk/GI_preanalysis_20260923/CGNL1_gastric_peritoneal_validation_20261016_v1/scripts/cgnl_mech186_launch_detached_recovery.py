from pathlib import Path
import subprocess,sys,json,datetime,hashlib
O=Path.cwd()/'mechanism_investigation_v2';root=O.parent;worker=root/'scripts/cgnl_mech185_resume_failed_files.py';wrapper=root/'scripts/cgnl_detached_wrapper_v1.py'
code='''from pathlib import Path
import json,sys,subprocess,datetime,traceback
p=Path(sys.argv[1]);m=json.loads(p.read_text(encoding='utf-8'));record=dict(m)
try:
 with Path(m['log']).open('wb') as f:
  proc=subprocess.Popen([m['python'],m['script']],cwd=m['cwd'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,close_fds=True);record['worker_pid']=proc.pid;Path(m['state']).write_text(json.dumps(record,indent=2),encoding='utf-8');record['exit_code']=proc.wait()
except Exception:
 record['launcher_error']=traceback.format_exc();record['exit_code']=1
record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();Path(m['receipt']).write_text(json.dumps(record,indent=2),encoding='utf-8')
'''
if wrapper.exists():assert wrapper.read_text(encoding='utf-8')==code
else:wrapper.write_text(code,encoding='utf-8')
assert not (O/'metadata/DETACHED_RECOVERY_186.json').exists();stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');stem=worker.stem+'_'+stamp;m={'python':sys.executable,'script':str(worker),'script_sha256':hashlib.sha256(worker.read_bytes()).hexdigest(),'cwd':str(root),'log':str(root/'logs'/(stem+'.log')),'state':str(root/'logs'/(stem+'.running.json')),'receipt':str(root/'logs'/(stem+'.exit.json')),'start_utc':stamp,'detach_reason':'old177 confirmed absent without completion after proxy andSSHconnection loss; onlystopped recovery resumed; otherliveproducers retained'};launch=root/'logs'/(stem+'.launch.json');launch.write_text(json.dumps(m,indent=2),encoding='utf-8')
proc=subprocess.Popen([sys.executable,str(wrapper),str(launch)],cwd=root,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True,close_fds=True);m.update(wrapper_pid=proc.pid,launch_record=str(launch));(O/'metadata/DETACHED_RECOVERY_186.json').write_text(json.dumps(m,indent=2),encoding='utf-8');print(json.dumps(m))
