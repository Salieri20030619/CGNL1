from pathlib import Path
import datetime
import fcntl
import hashlib
import json
import os
import select
import subprocess
import traceback

P = Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
F = P / 'continuous_exploration_20261023_v1/C04_full_fraction_qualification'
W = F / 'postqueue_annotation_20261025_v2'
Q = F / 'server_queue_20261025_v1'
W.mkdir(exist_ok=True)
lock = (W / 'runner.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

def save(name, value):
    value.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), whole_project_complete=False)
    (W / name).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')

try:
    assert not (W / 'EXIT.json').exists(), 'Prior terminal receipt exists; audit before rerun'
    script = P / 'scripts/cgnl_c04_annotate_joint_candidates.py'
    expected = '5b361d330db00121490c3139953ee9bc17f297da4464fd3d27b887dc4a5c8797'
    assert hashlib.sha256(script.read_bytes()).hexdigest() == expected
    try:
        fd = os.pidfd_open(2486092)
    except ProcessLookupError:
        assert (Q / 'EXIT.json').exists(), 'Queue absent without terminal receipt'
    else:
        assert 'cgnl_c04_server_queue_20261025.py' in Path('/proc/2486092/cmdline').read_bytes().decode()
        save('STATUS.json', dict(stage='waiting_native_queue_exit', pid=os.getpid(), queue_pid=2486092, annotator_sha256=expected))
        print('Waiting for original queue exit via pidfd', flush=True)
        select.select([fd], [], [])
        os.close(fd)
    state = json.loads((Q / 'STATUS.json').read_text())
    terminal = json.loads((Q / 'EXIT.json').read_text())
    joint = json.loads((Q / 'JOINT_RUN_EXIT.json').read_text())
    assert state['complete'] and len(state['files']) == 13
    assert all(x['stage'] in ['accepted', 'reused_complete', 'reused_pilot'] for x in state['files'].values())
    assert terminal['acquisition_search_complete'] and not terminal.get('error')
    assert joint['joint_exit'] == 0 and joint['R_exit'] == 0
    assert hashlib.sha256(script.read_bytes()).hexdigest() == expected
    output = F / 'joint_identification/CGNL1_SITE_ANNOTATIONS.json'
    assert not output.exists(), 'Existing annotations require review; do not overwrite'
    save('STATUS.json', dict(stage='annotating', pid=os.getpid(), new_disease_tests=0))
    result = subprocess.run(['python3', str(script)], capture_output=True, text=True)
    save('ANNOTATOR_EXIT.json', dict(exit=result.returncode, stdout=result.stdout, stderr=result.stderr))
    assert result.returncode == 0 and output.exists()
    parsed = json.loads(output.read_text())
    assert parsed['new_biological_tests'] == 0 and parsed['precursor_linkage_quantification_and_protein_gates_pending']
    save('EXIT.json', dict(annotation_complete=True, output=str(output), output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(), downstream_qualification_pending=True))
    save('STATUS.json', dict(stage='annotation_complete_downstream_pending', pid=os.getpid()))
    print('Annotation complete; linkage, quantification and protein gates remain pending', flush=True)
except Exception:
    save('FAILURE.json', dict(error=traceback.format_exc(), annotation_complete=False))
    raise
