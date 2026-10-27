
from pathlib import Path
import json,time,datetime,os
p=Path.cwd();m=p/'mechanism_investigation_v2'
names=['GSE286057_RAW_QA_LAUNCH_211.json','RESUME_PIPELINE_LAUNCH_192.json','GSE286057_RECOVERY_LAUNCH_209.json','GSE286057_ANALYSIS_LAUNCH_198c.json','DRA_QA_LAUNCH_200b.json','GSE286057_STATISTICS_LAUNCH_204c.json']
def read(f):
 try:return json.loads(f.read_text(encoding='utf-8'))
 except (FileNotFoundError,json.JSONDecodeError):return None
def save(f,v):
 tmp=f.with_name(f.name+'.writing');tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(f)
previous=None
while True:
 jobs=[]
 for name in names:
  launch=read(m/'metadata'/name)
  if not launch:jobs.append({'record':name,'state':'missing_launch'});continue
  receipt=read(Path(launch['receipt']));exists=Path('/proc',str(launch['wrapper_pid'])).exists()
  state='exited' if receipt is not None else ('running_or_waiting' if exists else 'missing_exit_receipt')
  jobs.append({'record':name,'state':state,'exit_code':receipt.get('exit_code') if receipt is not None else None,'process_exists':exists,'receipt':launch['receipt'],'log':launch['log']})
 quant=read(m/'metadata/DRA_QUANTIFICATION_STATUS.json') or {};raw=m/'raw/GSE286057_comparator_195'
 downloads=read(m/'metadata/GSE286057_DOWNLOAD_RECEIPTS_195.json') or []
 qa=read(m/'results/DRA_FULL_INDEPENDENT_QA_200.json') or {}
 snapshot={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'monitor_pid':os.getpid(),'interval_seconds':900,'quantified_libraries':len(quant.get('done',{})),'expected_libraries':quant.get('expected',12),'quantification_failures':quant.get('fail',{}),'DRA_independent_QA_pass':qa.get('pass'),'atlas_final_files':len(list(raw.glob('*.gz'))),'atlas_expected_files':108,'atlas_download_receipts':len(downloads),'atlas_failed_receipts':(0 if (read(m/'metadata/GSE286057_DOWNLOAD_RECOVERED_COMPLETE_209.json') or {}).get('all_pass') else sum(not r.get('success',False) for r in downloads)),'jobs':jobs,'project_goal_complete':False}
 snapshot['needs_review']=bool(any(j['state']=='missing_exit_receipt' or (j['state']=='exited' and j['exit_code']!=0) for j in jobs) or snapshot['atlas_failed_receipts'] or snapshot['quantification_failures'])
 save(m/'metadata/PROGRESS_LATEST_212.json',snapshot)
 with (m/'metadata/PROGRESS_HISTORY_212.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(snapshot,ensure_ascii=False)+chr(10))
 print(json.dumps(snapshot,ensure_ascii=False),flush=True)
 signature=[(j['record'],j['state'],j.get('exit_code')) for j in jobs]
 if previous is not None and signature!=previous:
  with (m/'metadata/PROGRESS_EVENTS_212.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(snapshot,ensure_ascii=False)+chr(10))
 previous=signature
 if all(j['state']=='exited' for j in jobs):
  save(m/'metadata/PROGRESS_MONITOR_COMPLETE_212.json',{'utc':snapshot['utc'],'all_tracked_jobs_terminal':True,'all_exit_codes_zero':all(j['exit_code']==0 for j in jobs),'output_integrity_and_scientific_interpretation_still_require_review':True,'project_goal_complete':False});break
 time.sleep(900)
