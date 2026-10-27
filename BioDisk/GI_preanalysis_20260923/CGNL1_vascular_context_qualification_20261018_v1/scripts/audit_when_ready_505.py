from pathlib import Path
import datetime,json,os,select,subprocess,traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,data):
    p=B/'metadata'/name;t=p.with_suffix(p.suffix+'.tmp')
    t.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8');t.replace(p)
with (B/'metadata/AUDIT505_LAUNCH.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),wrapper_pid=os.getpid(),started=now(),receipt=str(B/'metadata/AUDIT505.exit.json'),log=str(B/'logs/AUDIT505.log')),f,indent=2)
try:
    launch=json.loads((B/'metadata/CELL497_LAUNCH.json').read_text(encoding='utf-8'))
    receipt=B/'metadata/CELL497.exit.json'
    save('AUDIT505_STATE.json',dict(stage='native_completion_wait',dependency_pid=launch['pid']))
    if not receipt.exists():
        try:
            fd=os.pidfd_open(launch['pid']);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
        except ProcessLookupError:pass
    assert receipt.exists() and json.loads(receipt.read_text(encoding='utf-8'))['exit_code']==0,'Cell qualification failed or receipt missing'
    save('AUDIT505_STATE.json',dict(stage='numeric_audit',started=now()))
    result=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(B/'scripts/audit_pilot_500.py')])
    assert result.returncode==0,'Numeric audit failed'
    p=B/'results/pilot_numeric_audit_500/AUDIT500.json'
    assert json.loads(p.read_text(encoding='utf-8'))['numeric_checks_passed']
    status=dict(stage='ready_for_scientific_review',exit_code=0,finished=now(),result=str(p),mainline_validated=False)
except Exception:
    status=dict(stage='failed_needs_review',exit_code=1,finished=now(),error=traceback.format_exc())
save('AUDIT505_STATE.json',status)
save('AUDIT505.exit.json',status)

