from pathlib import Path
import os,json,datetime,select,time,errno
O=Path.cwd()/'mechanism_investigation_v2';det=json.loads((O/'metadata/DETACHED_RECOVERY_186.json').read_text(encoding='utf-8'));wanted=['cgnl_mech67_resumable_transfers.py','cgnl_mech46c_quantify_available_files.py','cgnl_mech91c_dra_target_when_ready.py','cgnl_mech94b_transcript_sensitivity_when_ready.py'];tracked=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  args=(p/'cmdline').read_bytes().split(b'\0');roles=[n for n in wanted if str(O.parent/'scripts'/n).encode() in args]
  if int(p.name)==det['wrapper_pid'] and det['launch_record'].encode() in args:roles.append('detached_recovery185_wrapper')
  if roles:tracked.append({'pid':int(p.name),'roles':roles})
 except (FileNotFoundError,PermissionError):pass
assert any('detached_recovery185_wrapper' in r['roles'] for r in tracked) or Path(det['receipt']).exists();fds={};native=True
try:
 for r in tracked:
  try:fds[os.pidfd_open(r['pid'])]=r
  except ProcessLookupError:pass
except (AttributeError,OSError) as e:
 if isinstance(e,OSError) and e.errno not in [errno.ENOSYS,errno.EINVAL]:raise
 native=False
 for f in fds:os.close(f)
 fds={}
state={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'tracked':tracked,'native_pidfd_wait':native,'read_only_process_observation':True,'remote_recovery_receipt':det['receipt']};(O/'metadata/PIPELINE_WATCH_187.json').write_text(json.dumps(state,indent=2),encoding='utf-8');s=json.loads((O/'STATUS.json').read_text(encoding='utf-8'));s.update(utc=state['utc'],last_goal_turn_classification='progress: SSHreconnectedandliveoriginalworkersconfirmed; stopped177recoveredasdetached185withverifiedchunkreuseanddirectroutefallback; nativeprocess-exitmonitorarmed');s['transport_recovery'].update(controller='cgnl_mech185_resume_failed_files.py',detached_launch='metadata/DETACHED_RECOVERY_186.json',old177_terminal_without_complete=True,new_failed_file_snapshot='metadata/LOCAL_FAILED_FILE_SNAPSHOT_186.json',owned_proxy_alias='sgs10-codex',proxy_restore_probe='metadata/ENA_PROXY_PROBE_183.json');s['pipeline_event_monitor']='metadata/PIPELINE_WATCH_187.json';(O/'STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8');print('Native process monitor armed',json.dumps(state),flush=True)
if native and fds:
 ready,_,_=select.select(list(fds),[],[]);events=[fds[f] for f in ready]
else:
 while True:
  events=[r for r in tracked if not (Path('/proc')/str(r['pid'])).exists()]
  if events or not tracked:break
  time.sleep(900)
for fd in fds:os.close(fd)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exited_processes':events,'process_exit_not_success_without_receipt':True,'markers':{}}
for p in [Path(det['receipt']),O/'metadata/FAILED_FILE_RECOVERY_COMPLETE_185.json',O/'metadata/DRA_QUANTIFICATION_COMPLETE.json',O/'results/DRA_TARGET_ANALYSIS_COMPLETE.json',O/'results/DRA_TRANSCRIPT_SENSITIVITY_COMPLETE.json']:
 if p.is_file():result['markers'][p.name]=json.loads(p.read_text(encoding='utf-8'))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');(O/'metadata'/('PIPELINE_EVENT_187_'+stamp+'.json')).write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
