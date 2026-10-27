from common import *
import subprocess,os,time,concurrent.futures,fcntl,signal,traceback
lock=(Q/'QUEUE.lock').open('w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:print('queue already running');sys.exit(2)
now=lambda:time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());state={'id':C['id'],'pid':os.getpid(),'started':now(),'status':'running','tasks':{t['id']:{'status':'pending','deps':t['deps']} for t in C['tasks']}};dump(Q/'STATE.json',state);dump(Q/'PID.json',{'pid':os.getpid(),'queue':str(Q)})
def execute(task):
 name=task['id'];D=Q/name;D.mkdir(exist_ok=True);env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2');start=now();code=1;pid=None;error=None
 try:
  with (D/'run.log').open('w') as log:
   p=subprocess.Popen([sys.executable,str(Q/'scripts'/task['script']),str(Q)],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env,start_new_session=True);pid=p.pid;dump(D/'STARTED.json',{'pid':pid,'started':start,'script':task['script']})
   try:code=p.wait(timeout=8*3600)
   except subprocess.TimeoutExpired:
    os.killpg(pid,signal.SIGTERM)
    try:p.wait(timeout=30)
    except subprocess.TimeoutExpired:os.killpg(pid,signal.SIGKILL);p.wait()
    error='task timed out';code=124
  if code==0 and not (D/'RESULT.json').exists():code=1;error='missing RESULT.json'
 except Exception:code=1;error=traceback.format_exc()
 result=json.loads((D/'RESULT.json').read_text()) if (D/'RESULT.json').exists() else {};status='failed' if code else ('skipped_condition' if result.get('status','').startswith(('skipped','insufficient')) else 'completed_pending_review');out={'status':status,'pid':pid,'exit_code':code,'started':start,'finished':now(),'result_status':result.get('status'),'error':error};dump(D/'EXIT.json',out)
 if code==0:dump(D/'COMPLETE.json',out)
 return out
terminal={'completed_pending_review','failed','skipped_condition','dependency_failed'};running={}
with concurrent.futures.ThreadPoolExecutor(max_workers=C['max_workers']) as pool:
 while True:
  for task in C['tasks']:
   name=task['id'];node=state['tasks'][name]
   if node['status']!='pending':continue
   deps=[state['tasks'][d]['status'] for d in task['deps']]
   if not all(s in terminal for s in deps):continue
   if any(s in ['failed','dependency_failed'] for s in deps) and not task.get('always_run'):
    node.update(status='dependency_failed',reason='upstream computational failure');continue
   if len(running)>=C['max_workers']:break
   node['status']='running';running[pool.submit(execute,task)]=name
  dump(Q/'STATE.json',state)
  if not running:
   if all(x['status'] in terminal for x in state['tasks'].values()):break
   raise RuntimeError('dependency deadlock')
  done,_=concurrent.futures.wait(running,return_when=concurrent.futures.FIRST_COMPLETED)
  for f in done:
   name=running.pop(f);state['tasks'][name].update(f.result());print(name,state['tasks'][name]['status'],flush=True)
  dump(Q/'STATE.json',state)
state['status']='completed_with_failures' if any(v['status'] in ['failed','dependency_failed'] for v in state['tasks'].values()) else 'completed_pending_review';state['finished']=now();dump(Q/'STATE.json',state);dump(Q/'QUEUE_EXIT.json',{'exit_code':1 if state['status']=='completed_with_failures' else 0,'state':state['status'],'finished':state['finished']});print(state['status'],flush=True)
