from pathlib import Path
import os,json,hashlib,subprocess,datetime,time
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
F=P/'continuous_exploration_20261023_v1/C04_full_fraction_qualification';W=F/'recovery_fraction11';W.mkdir(exist_ok=True)
name='phospho_fraction_11.raw';f=F/'raw'/(name+'.recovery.partial');prefix_bytes=187449344;prefix_sha='c70864b226624d8cbcbb28010e7bb66737c3a1356c66a6513551988709dba6aa'
def save(n,z):
    q=W/n;q.write_text(json.dumps(z,indent=2),encoding='utf-8');assert json.loads(q.read_text())==z
def digest(path,limit=None):
    h=hashlib.sha256();remaining=limit
    with path.open('rb') as handle:
        while remaining is None or remaining>0:
            b=handle.read(8*1024*1024 if remaining is None else min(remaining,8*1024*1024))
            if not b:break
            h.update(b)
            if remaining is not None:remaining-=len(b)
    assert remaining in (None,0)
    return h.hexdigest()
result={'pass':False,'stage':'starting','pid':os.getpid()}
try:
    lock=W/'ACTIVE.json';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
    save('ACTIVE.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':name})
    expected=next(x['bytes'] for x in json.loads((F/'ACQUISITION_SPEC_BEFORE_RESULTS.json').read_text())['files'] if x['name']==name)
    assert f.stat().st_size==prefix_bytes and digest(f)==prefix_sha,'Prefix transfer verification failed'
    assert not (F/'raw'/name).exists() and not (F/'raw'/(name+'.uploading')).exists()
    save('PREFIX_ACCEPTANCE.json',{'bytes':prefix_bytes,'sha256':prefix_sha,'matches_windows_prefix':True,'expected_complete_bytes':expected,'remaining_bytes':expected-prefix_bytes})
    print('Prefix verified; resuming remaining bytes with 7200-second deadline',flush=True)
    cmd=['curl','-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time','7200','--speed-limit','1024','--speed-time','180','--silent','--show-error','--continue-at','-','--output',str(f),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name]
    started=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True)
    save('TRANSFER_EXIT.json',{'exit':r.returncode,'stderr':r.stderr,'elapsed_seconds':time.monotonic()-started,'bytes':f.stat().st_size,'expected_bytes':expected,'certificate_validation_enabled':True})
    result['stage']='transfer_finished';assert r.returncode==0 and f.stat().st_size==expected,'Transfer incomplete; partial retained'
    assert digest(f,prefix_bytes)==prefix_sha
    with f.open('rb') as h:assert h.read(2)==b'\x01\xa1'
    sha=digest(f);save('RAW_ACQUISITION_RECEIPT.json',{'name':name,'bytes':expected,'sha256':sha,'prefix_verified':True,'transfer_exit':0})
    assert not (F/'raw'/name).exists() and not (F/'raw'/(name+'.uploading')).exists(),'Another writer appeared; preserve recovery file'
    f.rename(F/'raw'/(name+'.uploading'))
    with (W/'conversion_and_search.log').open('w',encoding='utf-8') as log:
        r=subprocess.run(['python3',str(P/'scripts/cgnl_c04_accept_convert_file.py'),name,sha],stdout=log,stderr=subprocess.STDOUT)
    result['stage']='conversion_and_search_finished';result['pipeline_exit']=r.returncode;assert r.returncode==0
    complete=F/'searches'/Path(name).stem/'FILE_SEARCH_COMPLETE.json';assert complete.exists()
    result.update(complete_receipt=str(complete),raw_sha256=sha)
    result['pass']=True
except Exception as e:
    result['error']=repr(e)
finally:
    result['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save('RECOVERY_EXIT.json',result);print(json.dumps(result),flush=True)
raise SystemExit(0 if result['pass'] else 1)
