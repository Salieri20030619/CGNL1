from pathlib import Path
import os,select,json,time,base64,hashlib,ctypes
def open_pidfd(pid):
 libc=ctypes.CDLL(None,use_errno=True)
 fd=libc.syscall(434,ctypes.c_int(pid),ctypes.c_uint(0))
 if fd<0:
  err=ctypes.get_errno()
  if err==3:raise ProcessLookupError(pid)
  raise OSError(err,os.strerror(err))
 return fd
R=Path(__file__).parent
jobs=[('wave01_proxy','wave01_core'),('wave02','wave02_extended'),('wave03','wave03_supplements')]
def save(n,x):
 p=R/n;t=R/(n+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def check(label,wave):
 try:
  ex=json.loads((R/(label+'.EXIT.json')).read_text());c=json.loads((R/wave/'COMPLETE.json').read_text());m=json.loads((R/wave/'MANIFEST.json').read_text())
  errs=[]
  for x in m:
   if x['status'] not in ['downloaded_verified','verified_existing']:errs.append(x['name']+': '+x.get('error','not verified'));continue
   p=Path(x['path'])
   if not p.exists() or p.stat().st_size!=x['bytes']:errs.append(x['name']+': file size changed')
   if p.suffix in ['.h5seurat','.h5ad','.h5']:
    try:
     import h5py
     with h5py.File(p,'r') as f:
      if not len(f.keys()):errs.append(x['name']+': HDF5 empty')
    except ImportError:pass
    except Exception as e:errs.append(x['name']+': HDF5 '+str(e))
   if x['dataset']=='GTEx_v8':
    refs=json.loads((R/'gtex_files.json').read_text()).get('items',[])
    z=next((y for y in refs if y['name'].endswith('/'+x['name'])),None)
    if z and z.get('md5Hash'):
     h=hashlib.md5()
     with p.open('rb') as f:
      for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
     if h.digest()!=base64.b64decode(z['md5Hash']):errs.append(x['name']+': publisher MD5 mismatch')
  return {'wave':wave,'process_exit':ex['exit_code'],'completion':c,'integrity_errors':errs,'status':'verified_complete' if ex['exit_code']==0 and c['exit_code']==0 and not errs else 'needs_attention'}
 except Exception as e:return {'wave':wave,'status':'needs_attention','error':str(e)}
poller=select.poll();fds={};results=[]
for label,wave in jobs:
 try:
  pid=int((R/(label+'.pid')).read_text());fd=open_pidfd(pid);poller.register(fd,select.POLLIN);fds[fd]=(label,wave)
 except ProcessLookupError:results.append(check(label,wave))
save('BATCH_COMPLETION_STATUS.json',{'status':'waiting_on_process_exit','results':results,'pending':[v[1] for v in fds.values()]})
while fds:
 for fd,event in poller.poll():
  label,wave=fds.pop(fd);poller.unregister(fd);os.close(fd);results.append(check(label,wave))
  save('BATCH_COMPLETION_STATUS.json',{'status':'waiting_on_process_exit' if fds else 'finished','results':results,'pending':[v[1] for v in fds.values()]})
save('DOWNLOAD_BATCHES_COMPLETE.json',{'exit_code':0 if all(x['status']=='verified_complete' for x in results) else 1,'results':results,'time':time.time(),'scope':'three queued batches only; ATAC raw and MR full-region resource gaps remain separately recorded'})
