import pathlib,sys,json,ctypes,os,select,subprocess,time
O=pathlib.Path(sys.argv[1]);jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text());libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(os.O_NONBLOCK);assert fd>=0;wd=libc.inotify_add_watch(fd,str(O/'raw/liver').encode(),0x00000008|0x00000080|0x00000100);assert wd>=0
last=-1
while True:
 ready=[]
 for j in jobs:
  p=O/'raw/liver'/(j['sample']+'.complete.json');data=O/'raw/liver'/(j['sample']+'.tar.gz')
  try:
   r=json.loads(p.read_text());ok=r.get('ok') and data.exists() and data.stat().st_size==r.get('bytes')
  except (FileNotFoundError,json.JSONDecodeError):ok=False
  if ok:ready.append(j['sample'])
 if len(ready)!=last:print('Verified input completion events:',len(ready),'/24',flush=True);last=len(ready)
 if len(ready)==len(jobs):break
 events=select.select([fd],[],[],900)[0]
 if events:
  try:os.read(fd,1024*1024)
  except BlockingIOError:pass
os.close(fd)
with (O/'results/full_screen.log').open('w') as log:
 r=subprocess.run([sys.executable,str(O/'scripts/screen_liver.py'),str(O)],stdout=log,stderr=subprocess.STDOUT)
out={'exit_code':r.returncode,'completion_marker':(O/'results/screen/COMPLETE.json').exists(),'log':str(O/'results/full_screen.log')};(O/'results/SCREEN_RUN_EXIT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True);sys.exit(r.returncode)
