from pathlib import Path
import ctypes,ctypes.wintypes as wt,json,os,sys,hashlib,datetime,math,subprocess
BASE=Path(r'C:\Users\31937\cgnl_transfer\restart_20261024')
D=BASE/'C04_full';W=D/'postqueue_recovery_v2';SOURCE=Path(r'C:\Users\31937\cgnl_c04_full_acquire.py')
ROOT='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1'
REMOTE=ROOT+'/continuous_exploration_20261023_v1/C04_full_fraction_qualification'
SSH=r'C:\Windows\System32\OpenSSH\ssh.exe';CURL=r'C:\Windows\System32\curl.exe';SCP=r'C:\Windows\System32\OpenSSH\scp.exe'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,obj):
    tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2),encoding='utf-8');tmp.replace(path)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def run(cmd,**kw):
    cmd=list(cmd)
    if cmd[0]=='scp':cmd[0]=SCP
    if cmd[0]=='ssh':cmd[0]=SSH
    return subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW,**kw)
def server_code(code):return run([SSH,'-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3','-'],input=code)
def deadline(remaining):return max(1800,min(28800,math.ceil(max(0,remaining)/8000*1.5)+300))

kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.OpenProcess.argtypes=[wt.DWORD,wt.BOOL,wt.DWORD];kernel.OpenProcess.restype=wt.HANDLE
kernel.WaitForSingleObject.argtypes=[wt.HANDLE,wt.DWORD];kernel.WaitForSingleObject.restype=wt.DWORD
kernel.CloseHandle.argtypes=[wt.HANDLE]
kernel.GetExitCodeProcess.argtypes=[wt.HANDLE,ctypes.POINTER(wt.DWORD)]
parent=22692;handle=None  # Original native session has already terminated
source=SOURCE.read_text(encoding='utf-8');marker='pending=[];errors=[]'
assert source.count(marker)==1
spec=json.loads((BASE/'C04_full_ACQUISITION_SPEC.json').read_text());assert len(spec['files'])==13
probe=json.loads((D/'CURRENT_CONNECTIVITY_PROBE.json').read_text());assert probe['pass']
assert (datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(probe['utc'])).total_seconds()<600
assert 'finished' in json.loads((D/'STATUS.json').read_text())
if '--check' in sys.argv:
    print(json.dumps({'parent_handle_live':bool(handle),'source_sha256':sha(SOURCE),'files':len(spec['files']),'one_additional_attempt_per_timeout':True,'sample_deadline_44MB':deadline(44317764),'max_deadline_seconds':28800,'no_existing_process_stopped':True}))
    if handle:kernel.CloseHandle(handle)
    raise SystemExit(0)
W.mkdir(exist_ok=True)
lock=W/'RUN_CHECKPOINT.json'
with lock.open('x',encoding='utf-8') as f:json.dump({'pid':os.getpid(),'parent_pid':parent,'utc':now(),'source_sha256':sha(SOURCE),'phase':'waiting_original_queue','native_wait':True},f,indent=2)
print('Original queue terminal verified; recovering failed inputs after fresh connectivity checks',flush=True)
if handle:
    assert kernel.WaitForSingleObject(handle,0xFFFFFFFF)==0
    ec=wt.DWORD();assert kernel.GetExitCodeProcess(handle,ctypes.byref(ec));kernel.CloseHandle(handle)
else:ec=wt.DWORD(0)
original=json.loads((D/'STATUS.json').read_text());assert 'finished' in original,'No terminal original queue receipt; inspect, do not restart'
save(W/'ORIGINAL_TERMINAL_STATUS.json',original)
save(W/'CONNECTIVITY_BEFORE_RESUME.json',probe)
ns={};exec(compile(source.split(marker)[0],str(SOURCE),'exec'),ns)
state={'utc_started':now(),'phase':'recovering','original_parent_exit':ec.value,'files':{},'complete':False}
def update(name,**kw):
    state['files'].setdefault(name,{}).update(kw);save(W/'STATUS.json',state);print(json.dumps({'file':name,**kw}),flush=True)
ns['update']=update;ns['run']=run
for item in ns['files']:
    name=item['name'];prior=original['files'].get(name,{});stage=prior.get('stage')
    if stage=='accepted':update(name,stage='reused_original_acceptance');continue
    if name=='phospho_fraction_11.raw':
        code=f"from pathlib import Path\nimport json\np=Path({REMOTE!r})\nf=p/'recovery_fraction11/RECOVERY_EXIT.json'\nz=json.loads(f.read_text()) if f.exists() else {{'pass':False,'pending':True}}\nz['search_complete']=(p/'searches/phospho_fraction_11/FILE_SEARCH_COMPLETE.json').exists()\nprint(json.dumps(z))\n"
        r=server_code(code)
        if r.returncode: update(name,stage='server_recovery_check_failed',stderr=r.stderr[-500:]);continue
        z=json.loads(r.stdout.strip().splitlines()[-1]);save(W/'FRACTION11_SERVER_STATE.json',z)
        update(name,stage='reused_server_acceptance' if z.get('pass') and z['search_complete'] else 'server_recovery_unresolved');continue
    if stage not in ['download_failed','transfer_failed']:
        update(name,stage='requires_failure_diagnosis',original_stage=stage);continue
    if stage=='download_failed' and not all(a.get('exit') in [6,28,56] for a in prior.get('attempts',[])):
        update(name,stage='non_transport_failure_requires_diagnosis');continue
    full=D/name;part=D/(name+'.partial');receipt=D/(name+'.receipt.json')
    try:
        if full.exists():
            rec=json.loads(receipt.read_text());assert full.stat().st_size==item['bytes'] and sha(full)==rec['sha256']
        else:
            before=part.stat().st_size if part.exists() else 0;assert before<=item['bytes'];seconds=deadline(item['bytes']-before)
            update(name,stage='resuming_once',prefix_bytes=before,remaining_bytes=item['bytes']-before,deadline_seconds=seconds)
            cmd=[CURL,'-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time',str(seconds),'--speed-limit','1024','--speed-time','180','--silent','--show-error']
            if before:cmd+=['--continue-at','-']
            cmd+=['--output',str(part),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name]
            r=run(cmd);attempt={'exit':r.returncode,'bytes_before':before,'bytes_after':part.stat().st_size if part.exists() else 0,'deadline_seconds':seconds,'stderr':r.stderr}
            save(W/(name+'.attempt.json'),attempt)
            if r.returncode!=0 or not part.exists() or part.stat().st_size!=item['bytes']:
                update(name,stage='bounded_recovery_incomplete',attempt=attempt)
                if r.returncode==6:
                    state['stop_reason']='DNS_failed_again; preserve remaining queue without repeated attempts';break
                continue
            with part.open('rb') as h:assert h.read(2)==b'\x01\xa1'
            digest=sha(part);part.rename(full);rec={'name':name,'bytes':item['bytes'],'expected':item['bytes'],'sha256':digest,'route':'native FTPS bounded postqueue resume','attempts':[attempt]};save(receipt,rec)
        if not ns['transfer_convert']((item,rec)):update(name,stage='recovery_pipeline_failed')
    except Exception as e:update(name,stage='recovery_exception',error=repr(e))
state['finished']=now();state['complete']=all(x.get('stage') in ['accepted','reused_original_acceptance','reused_server_acceptance'] for x in state['files'].values()) and len(state['files'])==12;save(W/'STATUS.json',state)
if state['complete']:
    code=f"from pathlib import Path\nimport subprocess,json\np=Path({ROOT!r})\nf=Path({REMOTE!r})\nr=subprocess.run(['python3',str(p/'scripts/cgnl_c04_joint_psm_q.py')],capture_output=True,text=True)\nz={{'joint_exit':r.returncode,'joint_stdout':r.stdout,'joint_stderr':r.stderr}}\nif r.returncode==0:\n q=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(p/'scripts/cgnl_c04_joint_q_independent.R'),str(f/'joint_identification')],capture_output=True,text=True);z.update(R_exit=q.returncode,R_stdout=q.stdout,R_stderr=q.stderr)\n(f/'POSTQUEUE_JOINT_RUN.json').write_text(json.dumps(z,indent=2),encoding='utf-8')\nprint(json.dumps(z))\n"
    r=server_code(code);save(W/'JOINT_TRIGGER.json',{'exit':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
print(json.dumps({'recovery_complete':state['complete'],'whole_project_complete':False}),flush=True)
