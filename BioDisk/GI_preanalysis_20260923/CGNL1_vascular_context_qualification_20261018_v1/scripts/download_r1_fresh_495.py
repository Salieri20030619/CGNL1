from pathlib import Path
import datetime, hashlib, json, os, subprocess, time, traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
manifest=json.loads((B/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text(encoding='utf-8'))
item=next(x for x in manifest['items'] if '_R1_' in x['name'])
D=Path(manifest['destination']).parent/'r1_fresh_495'
D.mkdir(exist_ok=False)
part=D/(item['name']+'.part')
part.touch(exist_ok=False)
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
state=dict(item,pid=os.getpid(),started=now(),route='server direct HTTPS',fresh_start=True,stage='downloading',qualified=False,attempts=[])
def save():
    p=B/'metadata/R1_FRESH495_STATE.json'
    t=p.with_suffix('.json.tmp')
    t.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8')
    t.replace(p)
save()
try:
    deadline=time.monotonic()+48*3600
    for attempt in range(1,121):
        if part.stat().st_size==item['bytes']: break
        if time.monotonic()>deadline: raise RuntimeError('48 hour retry limit reached')
        rec=dict(attempt=attempt,started=now(),start_bytes=part.stat().st_size)
        state['attempts'].append(rec)
        state['stage']='downloading'
        save()
        cmd=['curl','-q','--proxy','','--noproxy','*','--fail','--location','--connect-timeout','30','--max-time','14400','--speed-limit','1024','--speed-time','180','--retry','0','--continue-at','-','--output',str(part),item['url']]
        with (D/'curl.log').open('ab') as log:
            rec['curl_exit']=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=log).returncode
        rec.update(finished=now(),end_bytes=part.stat().st_size)
        save()
        if part.stat().st_size>item['bytes']: raise RuntimeError('Oversized download')
        if rec['curl_exit']==0 and part.stat().st_size==item['bytes']: break
        if rec['curl_exit'] in (23,26,33,60,63): raise RuntimeError('Non-transient curl error')
        state['stage']='retry_wait'
        save()
        time.sleep(60)
    state['actual_bytes']=part.stat().st_size
    if state['actual_bytes']!=item['bytes']: raise RuntimeError('Download incomplete')
    state['stage']='verifying'
    save()
    state['actual_md5']=subprocess.check_output(['md5sum',str(part)],text=True).split()[0]
    state['md5_matches']=state['actual_md5']==item['md5']
    g=subprocess.run(['gzip','-t',str(part)],capture_output=True,text=True)
    state.update(gzip_exit=g.returncode,gzip_stderr=g.stderr)
    state['qualified']=state['md5_matches'] and g.returncode==0
    if state['qualified']:
        target=D/item['name']
        assert not target.exists()
        part.rename(target)
        state['verified_path']=str(target)
except Exception:
    state['error']=traceback.format_exc()
state.update(stage='complete',finished=now())
save()

