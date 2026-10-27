from pathlib import Path
import datetime, json, os, select, subprocess, traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,data): (B/'metadata'/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
with (B/'metadata/CELL497_LAUNCH.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),started=now(),script=str(Path(__file__).resolve())),f,indent=2)
try:
    launch=json.loads((B/'metadata/PILOT496_LAUNCH.json').read_text(encoding='utf-8'))
    receipt=B/'metadata/PILOT496.exit.json'
    save('CELL497_STATE.json',dict(stage='waiting_for_quantification',dependency_pid=launch['pid']))
    if not receipt.exists():
        try:
            fd=os.pidfd_open(launch['pid'])
            poll=select.poll(); poll.register(fd,select.POLLIN); poll.poll(); os.close(fd)
        except ProcessLookupError: pass
    assert receipt.exists() and json.loads(receipt.read_text(encoding='utf-8'))['exit_code']==0, 'Quantification failed or exit receipt missing'
    source=(B/'scripts/qualify_pilot_cells_486.py').read_text(encoding='utf-8')
    source=source.replace('SPEC_485','SPEC_496').replace('PILOT485.exit','PILOT496.exit').replace('COMPLETE_485','COMPLETE_496')
    source=source.replace('pilot_cell_qualification_457','pilot_cell_qualification_497').replace('COMPLETE_457','COMPLETE_497')
    assert '485' not in source and '457' not in source
    output=B/'results/pilot_cell_qualification_497'
    assert not output.exists(), 'Preserve existing results'
    script=B/'scripts/qualify_pilot_cells_497.py'
    compile(source,str(script),'exec')
    with script.open('x',encoding='utf-8') as f: f.write(source)
    save('CELL497_STATE.json',dict(stage='qualifying_cells',started=now()))
    result=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(script)])
    assert result.returncode==0, 'Cell qualification failed'
    marker=output/'CELL_QUALIFICATION_COMPLETE_497.json'
    assert marker.is_file()
    save('CELL497.exit.json',dict(exit_code=0,finished=now(),result=str(marker)))
    save('CELL497_STATE.json',dict(stage='complete',finished=now(),result=str(marker)))
except Exception:
    error=traceback.format_exc()
    save('CELL497.exit.json',dict(exit_code=1,finished=now(),error=error))
    save('CELL497_STATE.json',dict(stage='failed',finished=now(),error=error))
    raise

