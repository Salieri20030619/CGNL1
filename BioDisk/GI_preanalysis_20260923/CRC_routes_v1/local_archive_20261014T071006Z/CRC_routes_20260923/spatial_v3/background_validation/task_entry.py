import pathlib,sys,json,subprocess,os,time,traceback,fcntl
R=pathlib.Path(sys.argv[1]);job=sys.argv[2];D=R/job;D.mkdir(exist_ok=True);lock=(D/'RUN.lock').open('w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:print('already running');sys.exit(2)
def write(name,j):
 p=D/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2));t.replace(p)
write('STARTED.json',{'job':job,'pid':os.getpid(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':'running'})
script={'V01':'v01_strict_masks.py','DATA_PREP':'data_prepare.py'}[job];env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
try:
 with (D/'run.log').open('w') as log:p=subprocess.run([sys.executable,str(R/'scripts'/script),str(R)],stdout=log,stderr=subprocess.STDOUT,env=env)
 code=p.returncode
except Exception:code=1;(D/'launcher_error.txt').write_text(traceback.format_exc())
result={'job':job,'pid':os.getpid(),'exit_code':code,'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':'completed_pending_review' if code==0 else 'failed'};write('EXIT.json',result)
if code==0:write('COMPLETE.json',result)
sys.exit(code)
