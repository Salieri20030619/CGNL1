from pathlib import Path
import subprocess,json,hashlib,time,datetime,shlex,ctypes
O=Path('D:/CGNL1_DRA_public_data_20261016');plan=json.loads((O/'PLAN.json').read_text(encoding='utf-8'));records={};attempts={}
R='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_gastric_peritoneal_validation_20261016_v1/mechanism_investigation_v2';py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
def save():
 p=O/'HANDOFF_STATUS.json';tmp=p.with_suffix('.writing');tmp.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'records':records,'attempts':attempts},indent=2),encoding='utf-8');tmp.replace(p)
def run(cmd,timeout=120,input=None):return subprocess.run(cmd,input=input,capture_output=True,timeout=timeout,creationflags=subprocess.CREATE_NO_WINDOW)
def remote(code,payload):
 r=run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10',py+' -c '+shlex.quote(code)],input=json.dumps(payload).encode('utf-8'));assert r.returncode==0,r.stderr.decode('utf-8');return json.loads(r.stdout)
def alive():
 k=ctypes.WinDLL('kernel32',use_last_error=True);k.OpenProcess.restype=ctypes.c_void_p;k.WaitForSingleObject.argtypes=[ctypes.c_void_p,ctypes.c_ulong];k.CloseHandle.argtypes=[ctypes.c_void_p];h=k.OpenProcess(0x100000,False,plan['pid'])
 if not h:raise OSError(ctypes.get_last_error(),'Cannot observe download process; not proof of exit')
 result=k.WaitForSingleObject(h,0);k.CloseHandle(h);assert result in [0,258];return result==258
prepare="""from pathlib import Path
import json,sys,hashlib
d=json.load(sys.stdin);O=Path(d['root']);p=O/'raw/local_verified_D';p.mkdir(exist_ok=True);f=p/d['name'];valid=False
if f.exists():
 h=hashlib.md5()
 with f.open('rb') as s:
  while b:=s.read(4194304):h.update(b)
 assert f.stat().st_size==d['bytes'] and h.hexdigest()==d['md5'];valid=True
print(json.dumps({'stage':str(f),'already_verified':valid}))
"""
verify="""from pathlib import Path
import json,sys,hashlib,os,datetime
d=json.load(sys.stdin);O=Path(d['root']);stage=O/'raw/local_verified_D'/d['name'];p=stage if stage.exists() else stage.with_name(stage.name+'.incoming');h=hashlib.md5()
with p.open('rb') as s:
 while b:=s.read(4194304):h.update(b)
assert p.stat().st_size==d['bytes'] and h.hexdigest()==d['md5']
if p!=stage:os.link(p,stage)
dst=O/'raw/DRA014736'/d['name'];created=False
try:os.link(stage,dst);created=True
except FileExistsError:
 h2=hashlib.md5()
 with dst.open('rb') as s:
  while b:=s.read(4194304):h2.update(b)
 assert dst.stat().st_size==d['bytes'] and h2.hexdigest()==d['md5']
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'name':d['name'],'bytes':d['bytes'],'md5':h.hexdigest(),'master_hardlink_created':created,'controller_receipts_unchanged':True,'success':True};(O/'metadata'/('LOCAL_IMPORT_'+d['name']+'.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt))
"""
while True:
 rp=O/'RECEIPTS.json';receipts=json.loads(rp.read_text(encoding='utf-8')) if rp.exists() else []
 for d in receipts:
  name=d['name']
  if not d.get('success') or records.get(name,{}).get('success') or attempts.get(name,0)>=3:continue
  attempts[name]=attempts.get(name,0)+1
  try:
   p=O/name;h=hashlib.md5()
   with p.open('rb') as s:
    while b:=s.read(4194304):h.update(b)
   assert p.stat().st_size==d['bytes'] and h.hexdigest()==d['md5'];payload=dict(d,root=R);stage=remote(prepare,payload)
   if not stage['already_verified']:
    r=run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20',str(p),'sgs10:'+stage['stage']+'.incoming'],timeout=10800);assert r.returncode==0,r.stderr.decode('utf-8')
   records[name]=remote(verify,payload)
  except Exception as e:records[name]={'success':False,'error':repr(e),'attempt':attempts[name]}
  save();print(name,records[name],flush=True)
 done=O/'COMPLETE.json'
 if done.exists() and all(not d.get('success') or records.get(d['name'],{}).get('success') or attempts.get(d['name'],0)>=3 for d in receipts):break
 if not done.exists() and not alive():raise RuntimeError('Download exited without completion receipt')
 time.sleep(900)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':len(records)==22 and all(r.get('success') for r in records.values()),'records':records};(O/'HANDOFF_COMPLETE.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
