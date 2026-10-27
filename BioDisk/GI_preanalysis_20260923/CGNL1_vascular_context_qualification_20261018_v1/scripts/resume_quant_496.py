from pathlib import Path
import os, json, gzip, collections, time, subprocess, datetime, traceback
b=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,data): (b/'metadata'/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
with (b/'metadata/PILOT496_LAUNCH.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),started=now(),script=str(Path(__file__).resolve())),f,indent=2)
try:
    r1=json.loads((b/'metadata/R1_FRESH495_STATE.json').read_text(encoding='utf-8'))
    prior=json.loads((b/'metadata/DIRECT494_STATE.json').read_text(encoding='utf-8'))
    r2=next(x for x in prior['results'] if '_R2_' in x['name'])
    manifest=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text(encoding='utf-8'))
    files={}
    for receipt in (r1,r2):
        expected=next(x for x in manifest['items'] if x['name']==receipt['name'])
        assert receipt['qualified'] and receipt['md5_matches'] and receipt['gzip_exit']==0
        assert receipt['actual_md5']==expected['md5']
        p=Path(receipt['verified_path'])
        assert p.is_file() and p.stat().st_size==expected['bytes']
        files['R1' if '_R1_' in p.name else 'R2']=p
    save('PILOT_DOWNLOAD_COMPLETE_496.json',dict(utc=now(),files=[r1,r2],all_md5_verified=True))
    spec=json.loads((b/'metadata/PILOT_QUANTIFICATION_SPEC_485.json').read_text(encoding='utf-8'))
    spec['output']=str(Path(spec['output']).with_name('pilot_quantification_496'))
    spec['download_dependency']=str(b/'metadata/PILOT_DOWNLOAD_COMPLETE_496.json')
    assert not Path(spec['output']).exists(), 'Output exists; do not overwrite'
    save('PILOT_QUANTIFICATION_SPEC_496.json',spec)
    source=(b/'scripts/validate_and_quantify_pilot_485.py').read_text(encoding='utf-8')
    assert source.count('counts=[collections.Counter()')==1
    body=source[source.index('counts=[collections.Counter()'):].replace('_485','_496')
    save('PILOT496.running.json',dict(pid=os.getpid(),stage='paired_FASTQ_QA_then_alevin',started=now()))
    exec(compile(body,'validated_pipeline_485_body_version496','exec'))
    save('PILOT496.exit.json',dict(exit_code=0,finished=now()))
except Exception:
    save('PILOT496.exit.json',dict(exit_code=1,finished=now(),error=traceback.format_exc()))
    raise

