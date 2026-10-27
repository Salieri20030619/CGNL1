from pathlib import Path
import json, subprocess, concurrent.futures, datetime, traceback, os

b = Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
m = json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text(encoding='utf-8'))
root = Path(m['destination']).parent
stage = root/'integrity_repair_491'
exe = json.loads((b/'metadata/ASPERA_CLIENT_475.json').read_text(encoding='utf-8'))['executable']
key = root/'software/aspera_474/ibm_public_client_rsa.pem'
def save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()
def run(item):
    target = stage/item['name']
    logs = stage/('logs_'+item['name'])
    logs.mkdir()
    receipt = dict(item, started=now(), qualified=False)
    rp = b/'metadata'/('REPAIR491_'+item['name']+'.json')
    try:
        subprocess.run(['cp','--reflink=auto','--preserve=timestamps','--',str(root/'aspera_probe_480'/item['name']),str(target)], check=True)
        receipt['stage']='full_source_checksum_compare_and_repair'
        save(rp, receipt)
        cmd=[exe,'-Q','-k','3','--overwrite=diff','--file-checksum=md5','--file-manifest=text','--file-manifest-path='+str(logs),'-l','20M','-P','33001','-i',str(key),'-L',str(logs),'era-fasp@fasp.sra.ebi.ac.uk:/vol1/run/ERR100/ERR10026294/'+item['name'],str(stage)]
        with (logs/'command.log').open('wb') as f:
            receipt['ascp_exit']=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=86400).returncode
        receipt['actual_bytes']=target.stat().st_size
        receipt['actual_md5']=subprocess.check_output(['md5sum',str(target)],text=True).split()[0]
        g=subprocess.run(['gzip','-t',str(target)],capture_output=True,text=True)
        receipt['gzip_exit']=g.returncode
        receipt['gzip_stderr']=g.stderr
        receipt['qualified']=receipt['ascp_exit']==0 and receipt['actual_bytes']==item['bytes'] and receipt['actual_md5']==item['md5'] and g.returncode==0
    except Exception:
        receipt['error']=traceback.format_exc()
    receipt['finished']=now()
    save(rp,receipt)
    return receipt
stage.mkdir(exist_ok=False)
save(b/'metadata/REPAIR491_LAUNCH.json',dict(pid=os.getpid(),started=now(),script=str(Path(__file__).resolve())))
save(b/'metadata/REPAIR491_STATE.json',dict(stage='running',started=now(),original_files_preserved=True))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(run,m['items']))
save(b/'metadata/REPAIR491_STATE.json',dict(stage='complete',finished=now(),results=results,all_qualified=all(x['qualified'] for x in results)))
