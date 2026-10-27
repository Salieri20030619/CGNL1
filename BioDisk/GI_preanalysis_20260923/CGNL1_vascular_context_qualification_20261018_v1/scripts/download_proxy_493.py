from pathlib import Path
import concurrent.futures, datetime, fcntl, json, os, subprocess, time, traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
M=json.loads((B/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text(encoding='utf-8'))
D=Path(M['destination']).parent/'proxy_redownload_492'
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,data):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)
def download(item):
    result=dict(item,started=now(),qualified=False,stage='downloading',attempts=[])
    receipt=B/'metadata'/('PROXY493_'+item['name']+'.json')
    part=D/(item['name']+'.part')
    target=D/item['name']
    try:
        assert not target.exists(), 'Verified target already exists; inspect before resuming'
        assert part.exists() and part.stat().st_size<=item['bytes']
        result['resume_bytes']=part.stat().st_size
        save(receipt,result)
        deadline=time.monotonic()+48*3600
        for attempt in range(1,121):
            if part.stat().st_size==item['bytes']: break
            if time.monotonic()>deadline: raise RuntimeError('48 hour retry limit reached')
            command=['curl','--proxy','http://127.0.0.1:7897','--noproxy','','--fail','--location','--connect-timeout','30','--max-time','3600','--speed-limit','1024','--speed-time','180','--retry','0','--continue-at','-','--output',str(part),item['url']]
            record=dict(attempt=attempt,started=now(),start_bytes=part.stat().st_size)
            result['attempts'].append(record)
            result['stage']='downloading'
            save(receipt,result)
            with (D/(item['name']+'.resume493.log')).open('ab') as log:
                record['curl_exit']=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=log,stderr=log).returncode
            record.update(finished=now(),end_bytes=part.stat().st_size)
            save(receipt,result)
            if part.stat().st_size>item['bytes']: raise RuntimeError('Oversized partial file')
            if record['curl_exit']==0 and part.stat().st_size==item['bytes']: break
            if record['curl_exit'] in (23,26,33,60,63): raise RuntimeError('Non-transient curl error; retained partial')
            result['stage']='retry_wait'
            save(receipt,result)
            time.sleep(60)
        result['actual_bytes']=part.stat().st_size
        if result['actual_bytes']!=item['bytes']: raise RuntimeError('Incomplete after retry budget')
        result['stage']='verifying'
        save(receipt,result)
        result['actual_md5']=subprocess.check_output(['md5sum',str(part)],text=True).split()[0]
        result['md5_matches']=result['actual_md5']==item['md5']
        check=subprocess.run(['gzip','-t',str(part)],capture_output=True,text=True)
        result.update(gzip_exit=check.returncode,gzip_stderr=check.stderr)
        result['qualified']=result['md5_matches'] and check.returncode==0
        if result['qualified']:
            assert not target.exists()
            part.rename(target)
            result['verified_path']=str(target)
    except Exception:
        result['error']=traceback.format_exc()
    result.update(stage='complete',finished=now())
    save(receipt,result)
    return result
assert D.is_dir()
lock=(D/'resume493.lock').open('a')
fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (B/'metadata/PROXY493_LAUNCH.json').exists(), 'Launch already recorded; inspect state'
save(B/'metadata/PROXY493_LAUNCH.json',dict(pid=os.getpid(),started=now(),directory=str(D),proxy='http://127.0.0.1:7897',retry_limit=120,deadline_hours=48))
save(B/'metadata/PROXY493_STATE.json',dict(stage='running',started=now()))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(download,M['items']))
save(B/'metadata/PROXY493_STATE.json',dict(stage='complete',finished=now(),all_qualified=all(x['qualified'] for x in results),results=results))
