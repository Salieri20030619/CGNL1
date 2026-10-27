from pathlib import Path
import concurrent.futures, datetime, json, os, subprocess, traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
M=json.loads((B/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text(encoding='utf-8'))
D=Path(M['destination']).parent/'proxy_redownload_492'
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,data):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)
def download(item):
    result=dict(item,started=now(),qualified=False,stage='downloading')
    receipt=B/'metadata'/('PROXY492_'+item['name']+'.json')
    part=D/(item['name']+'.part')
    save(receipt,result)
    try:
        command=['curl','--proxy','http://127.0.0.1:7897','--noproxy','','--fail','--location','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--retry','10','--retry-delay','15','--continue-at','-','--output',str(part),item['url']]
        with (D/(item['name']+'.curl.log')).open('wb') as log:
            result['curl_exit']=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=log,stderr=log).returncode
        result['stage']='verifying'
        save(receipt,result)
        result['actual_bytes']=part.stat().st_size
        if result['curl_exit']!=0 or result['actual_bytes']!=item['bytes']:
            raise RuntimeError('Transfer failed or size mismatch; partial retained')
        result['actual_md5']=subprocess.check_output(['md5sum',str(part)],text=True).split()[0]
        result['md5_matches']=result['actual_md5']==item['md5']
        check=subprocess.run(['gzip','-t',str(part)],capture_output=True,text=True)
        result['gzip_exit']=check.returncode
        result['gzip_stderr']=check.stderr
        result['qualified']=result['md5_matches'] and check.returncode==0
        if result['qualified']:
            target=D/item['name']
            assert not target.exists()
            part.rename(target)
            result['verified_path']=str(target)
    except Exception:
        result['error']=traceback.format_exc()
    result['stage']='complete'
    result['finished']=now()
    save(receipt,result)
    return result
D.mkdir(exist_ok=False)
save(B/'metadata/PROXY492_LAUNCH.json',dict(pid=os.getpid(),started=now(),directory=str(D),proxy='http://127.0.0.1:7897',route='sgs10-codex reverse tunnel to Windows proxy; no local raw-data files'))
save(B/'metadata/PROXY492_STATE.json',dict(stage='running',started=now()))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(download,M['items']))
save(B/'metadata/PROXY492_STATE.json',dict(stage='complete',finished=now(),all_qualified=all(x['qualified'] for x in results),results=results))
