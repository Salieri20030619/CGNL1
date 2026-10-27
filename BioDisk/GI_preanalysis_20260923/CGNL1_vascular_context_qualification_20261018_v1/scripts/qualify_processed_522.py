from pathlib import Path
import collections,datetime,hashlib,json,os,subprocess,time,traceback
import h5py,numpy as np
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/independent_processed_522');D.mkdir(exist_ok=False)
manifest=json.loads((B/'metadata/external_sources_521/mendeley_files.json').read_text(encoding='utf-8'))
item=next(x for x in manifest if x['filename'].endswith('.h5ad'));part=D/(item['filename']+'.part')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
state=dict(pid=os.getpid(),started=now(),stage='downloading',dataset='10.17632/yd7bb3shn5.1',scope='Source qualification only; no HGP assignment inferred from expression',file_manifest=item)
def save():
    p=B/'metadata/PROCESSED522_STATE.json';t=p.with_suffix('.json.tmp');t.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8');t.replace(p)
save()
def decode(v):
    if isinstance(v,bytes):return v.decode('utf-8')
    if isinstance(v,np.generic):return v.item()
    return v
try:
    for attempt in range(1,4):
        with (D/'download.log').open('ab') as log:
            r=subprocess.run(['curl','-q','--proxy','','--noproxy','*','--fail','--location','--continue-at','-','--connect-timeout','30','--max-time','7200','--speed-time','180','--speed-limit','1024','--output',str(part),item['content_details']['download_url']],stdout=log,stderr=log)
        state.update(download_attempt=attempt,curl_exit=r.returncode);save()
        if r.returncode==0:break
        time.sleep(30)
    assert r.returncode==0 and part.stat().st_size==item['size'],'Transfer incomplete'
    h=hashlib.sha256()
    with part.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    state['actual_sha256']=h.hexdigest();assert h.hexdigest()==item['content_details']['sha256_hash'],'Source checksum mismatch'
    target=D/item['filename'];part.rename(target);state.update(stage='inspecting_metadata',verified_path=str(target));save()
    with h5py.File(target,'r') as f:
        state['root_keys']=list(f.keys());summaries={}
        for key,obj in f['obs'].items():
            if isinstance(obj,h5py.Group) and 'categories' in obj and 'codes' in obj:
                categories=[decode(v) for v in obj['categories'][()]];counts=collections.Counter(obj['codes'][()].tolist())
                summaries[key]=dict(encoding='categorical',counts={str(categories[i]) if i>=0 else '<missing>':int(n) for i,n in counts.items()})
            elif isinstance(obj,h5py.Dataset):
                values=obj[()]
                if values.dtype.kind in 'OSU':
                    counts=collections.Counter(decode(v) for v in values)
                    summaries[key]=dict(unique_values=len(counts),counts=dict(counts) if len(counts)<=100 else None,length=len(values))
                elif values.dtype.kind in 'biuf':
                    summaries[key]=dict(dtype=str(values.dtype),length=len(values),min=float(np.nanmin(values)),median=float(np.nanmedian(values)),max=float(np.nanmax(values)))
        state['obs_fields']=summaries;presence={}
        for key,obj in f['var'].items():
            if isinstance(obj,h5py.Dataset) and obj.dtype.kind in 'OSU':
                values=[decode(v) for v in obj[()]]
                hits=[i for i,v in enumerate(values) if v=='CGNL1']
                presence[key]=dict(entries=len(values),CGNL1_exact_matches=len(hits),match_indices=hits)
        state['var_target_presence']=presence
        x=f['X'];state['X']=dict(type=type(x).__name__,encoding=decode(x.attrs.get('encoding-type','')),shape=[int(v) for v in x.attrs.get('shape',x.shape if isinstance(x,h5py.Dataset) else [])],dtype=str(x.dtype) if isinstance(x,h5py.Dataset) else str(x['data'].dtype),stored_values=int(x['data'].size) if isinstance(x,h5py.Group) else int(x.size))
        state['layers']=list(f['layers'].keys()) if 'layers' in f else []
        state['uns_keys']=list(f['uns'].keys()) if 'uns' in f else []
    state.update(stage='qualified_file_pending_scientific_review',finished=now(),mainline_validated=False)
except Exception:state.update(stage='failed',error=traceback.format_exc(),finished=now())
save();print(json.dumps(state),flush=True)

