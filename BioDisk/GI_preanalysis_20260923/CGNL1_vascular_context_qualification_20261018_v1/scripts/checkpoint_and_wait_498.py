from pathlib import Path
import datetime, json, os, select, sys
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp')
    t.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');t.replace(p)
pointer=B.parent/'CGNL1_CURRENT_STATE_20261017.json'
backup=B/'metadata/GLOBAL_POINTER_BEFORE_498.json'
if not backup.exists():
    old=json.loads(pointer.read_text(encoding='utf-8'))
    with backup.open('x',encoding='utf-8') as f: json.dump(old,f,indent=2)
    checkpoint=dict(utc=now(),whole_project_complete=False,mainline_validated=False,
        stage='paired_FASTQ_QA_then_quantification_then_cell_qualification',
        download_evidence=[str(B/'metadata/R1_FRESH495_STATE.json'),str(B/'metadata/DIRECT494_STATE.json')],
        valid_R1_source='R1_FRESH495 only; earlier R1 copies failed integrity',
        valid_R2_source='DIRECT494 R2 result',
        quantification_launch=str(B/'metadata/PILOT496_LAUNCH.json'),
        cell_qualification_launch=str(B/'metadata/CELL497_LAUNCH.json'),
        next_review='Inspect exit receipts, output integrity, mapping and cell identity; a single donor supports descriptive feasibility only',
        automatic_whole_cohort_launch=False,prior_checkpoint=old.get('latest_checkpoint'))
    cp=B/'metadata/PROJECT_CHECKPOINT_498.json';save(cp,checkpoint)
    old.update(utc=now(),stage=checkpoint['stage'],latest_checkpoint=str(cp),
        pending_launches=[checkpoint['quantification_launch'],checkpoint['cell_qualification_launch']],
        next_native_completion=str(B/'metadata/CELL497.exit.json'),whole_project_complete=False,
        runtime_goal_state='active',runtime_goal_state_verified_utc=now(),
        qualified_download_evidence=checkpoint['download_evidence'])
    save(pointer,old)
if '--wait' in sys.argv:
    receipt=B/'metadata/CELL497.exit.json'
    launch=json.loads((B/'metadata/CELL497_LAUNCH.json').read_text(encoding='utf-8'))
    if not receipt.exists():
        try:
            fd=os.pidfd_open(launch['pid']); poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
        except ProcessLookupError: pass
    assert receipt.exists(), 'Process exited without receipt'
    print(receipt.read_text(encoding='utf-8'),flush=True)
else:
    print('Checkpoint 498 saved',flush=True)

