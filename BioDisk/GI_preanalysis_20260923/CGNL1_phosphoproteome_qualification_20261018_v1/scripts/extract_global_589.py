from pathlib import Path
import ctypes,ctypes.util,json,hashlib,traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H02')
dest=D/'global_589';dest.mkdir(exist_ok=True)
lib=ctypes.CDLL(ctypes.util.find_library('archive'))
def bind(name,args,restype):
 f=getattr(lib,name);f.argtypes=args;f.restype=restype;return f
ptr=ctypes.c_void_p
new=bind('archive_read_new',[],ptr)
filters=bind('archive_read_support_filter_all',[ptr],ctypes.c_int)
formats=bind('archive_read_support_format_all',[ptr],ctypes.c_int)
op=bind('archive_read_open_filename',[ptr,ctypes.c_char_p,ctypes.c_size_t],ctypes.c_int)
nex=bind('archive_read_next_header',[ptr,ctypes.POINTER(ptr)],ctypes.c_int)
name=bind('archive_entry_pathname',[ptr],ctypes.c_char_p)
size=bind('archive_entry_size',[ptr],ctypes.c_longlong)
kind=bind('archive_entry_filetype',[ptr],ctypes.c_uint)
data=bind('archive_read_data',[ptr,ptr,ctypes.c_size_t],ctypes.c_ssize_t)
error=bind('archive_error_string',[ptr],ctypes.c_char_p)
free=bind('archive_read_free',[ptr],ctypes.c_int)
a=new();filters(a);formats(a);records=[]
try:
 assert op(a,str(D/'MaxQuant_results_and_its_relative_supplementary_materials.rar').encode(),65536)==0
 e=ptr();buffer=ctypes.create_string_buffer(1024*1024);total=0
 while True:
  ret=nex(a,ctypes.byref(e))
  if ret==1:break
  assert ret==0,(ret,error(a))
  n=name(e).decode('utf-8');sz=size(e);typ=kind(e)
  if typ==0o040000:continue
  assert typ==0o100000,(n,typ)
  total+=sz;assert sz<4*1024**3 and total<12*1024**3
  p=(dest/n.replace('\\','/')).resolve();assert dest.resolve() in p.parents
  keep=p.suffix.lower() in ['.txt','.csv','.tsv','.xlsx','.xls','.xml','.pdf','.docx','.doc']
  f=None
  if keep:p.parent.mkdir(parents=True,exist_ok=True);f=p.open('wb')
  h=hashlib.sha256();count=0
  try:
   while True:
    got=data(a,buffer,len(buffer))
    assert got>=0,(n,got,error(a))
    if got==0:break
    chunk=buffer.raw[:got];h.update(chunk);count+=got
    if f:f.write(chunk)
  finally:
   if f:f.close()
  assert count==sz
  records.append(dict(member=n,bytes=count,sha256=h.hexdigest(),saved=str(p) if keep else None))
  print(n,count,flush=True)
 result=dict(exit_code=0,all_entries_read_without_libarchive_error=True,records=records)
except Exception:
 result=dict(exit_code=1,error=traceback.format_exc(),records=records);raise
finally:
 free(a)
 (B/'metadata/H02_GLOBAL_EXTRACTION589.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
