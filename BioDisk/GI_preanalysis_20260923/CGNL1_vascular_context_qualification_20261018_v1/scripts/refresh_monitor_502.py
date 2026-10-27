from pathlib import Path
import datetime,json
BASE=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
B=BASE/'CGNL1_vascular_context_qualification_20261018_v1'
registry=BASE/'CGNL1_pharmacology_20261017_v1/metadata/PROJECT_JOB_REGISTRY_216.json'
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');t.replace(p)
old=json.loads(registry.read_text(encoding='utf-8'))
with (B/'metadata/MONITOR_REGISTRY_BEFORE_502.json').open('x',encoding='utf-8') as f: json.dump(old,f,indent=2)
new_records=[]
for name in ('PILOT496','CELL497'):
    launch=json.loads((B/'metadata'/(name+'_LAUNCH.json')).read_text(encoding='utf-8'))
    adapter=dict(wrapper_pid=launch['pid'],receipt=str(B/'metadata'/(name+'.exit.json')),log=str(B/'logs'/(name+'.log')),original_launch=str(B/'metadata'/(name+'_LAUNCH.json')),pid_semantics='actual root worker process; no separate wrapper')
    assert Path('/proc',str(launch['pid'])).exists() or Path(adapter['receipt']).exists()
    path=B/'metadata'/(name+'_MONITOR_ADAPTER_502.json');save(path,adapter);new_records.append(str(path))
replaced=[p for p in old['jobs'] if Path(p).name in ('ASPERA484_LAUNCH.json','PILOT485_LAUNCH.json','CELL486_LAUNCH.json')]
old['jobs']=[p for p in old['jobs'] if p not in replaced]+new_records
old.setdefault('superseded_jobs',[]).extend(p for p in replaced if p not in old['superseded_jobs'])
old.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scientific_status='R1 fresh495 and R2 direct494 passed MD5/gzip. PILOT496 paired QA then quantification; CELL497 target-independent qualification. One-donor feasibility only; mainline unvalidated.',project_goal_complete=False)
for record in old['jobs']:
    adapter=json.loads(Path(record).read_text(encoding='utf-8'))
    assert all(k in adapter for k in ('wrapper_pid','receipt','log'))
save(registry,old)
print(json.dumps(dict(updated=True,new_jobs=new_records,superseded=replaced)))

