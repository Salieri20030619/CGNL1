from pathlib import Path
import json,datetime,time
O=Path.cwd()/'mechanism_investigation_v2';out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'processes':[],'completed_runs':[],'recovery_receipts':[],'markers':{},'active_ranges':[]}
watch=json.loads((O/'metadata/PIPELINE_WATCH_187.json').read_text(encoding='utf-8'));det=json.loads((O/'metadata/DETACHED_RECOVERY_186.json').read_text(encoding='utf-8'));state=json.loads(Path(det['state']).read_text(encoding='utf-8'));pids=[r['pid'] for r in watch['tracked']]+[state['worker_pid']]
for pid in pids:
 p=Path('/proc')/str(pid)
 if p.exists():
  try:out['processes'].append({'pid':pid,'state':next(l for l in (p/'status').read_text(encoding='utf-8').splitlines() if l.startswith('State:'))})
  except FileNotFoundError:out['processes'].append({'pid':pid,'state':'disappeared_during_read'})
 else:out['processes'].append({'pid':pid,'state':'absent'})
p=O/'metadata/DRA_QUANTIFICATION_STATUS.json'
if p.exists():
 q=json.loads(p.read_text(encoding='utf-8'));out['completed_runs']=list(q['done']);out['failed_runs']=q['fail']
p=O/'metadata/FAILED_FILE_RECOVERY_RECEIPTS_185.json'
if p.exists():out['recovery_receipts']=[{k:r.get(k) for k in ['name','success','end_utc','error']} for r in json.loads(p.read_text(encoding='utf-8'))]
for p in [Path(det['receipt']),O/'metadata/FAILED_FILE_RECOVERY_COMPLETE_185.json',O/'metadata/DRA_QUANTIFICATION_COMPLETE.json',O/'results/DRA_TARGET_ANALYSIS_COMPLETE.json']:
 if p.exists():out['markers'][p.name]=json.loads(p.read_text(encoding='utf-8'))
for folder in (O/'raw/recovery_resumed_185').glob('*.ranges'):
 bins=list(folder.glob('*.bin'));parts=list(folder.glob('*.partial'))
 if parts:out['active_ranges'].append({'file':folder.name,'verified_bytes':sum(p.stat().st_size for p in bins),'partial_bytes':sum(p.stat().st_size for p in parts),'newest_partial_age_seconds':min(time.time()-p.stat().st_mtime for p in parts)})
(O/'metadata/PIPELINE_INTERVAL_CHECK_189.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
