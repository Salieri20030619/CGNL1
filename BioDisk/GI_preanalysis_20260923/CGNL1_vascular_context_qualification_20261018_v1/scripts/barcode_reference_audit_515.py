from pathlib import Path
import collections,datetime,hashlib,json,os,subprocess,time,traceback,gzip,csv
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/technical_probe_513')
Q=D.parent/'pilot_quantification_496'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
state=dict(pid=os.getpid(),started=now(),stage='reference_download',biological_inference=False)
def save():
    p=B/'metadata/BARCODE_AUDIT_515.json';t=p.with_suffix('.json.tmp');t.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8');t.replace(p)
save()
try:
    p=D/'3M-february-2018.malva.txt';url='https://bimsbstatic.mdc-berlin.de/rajewsky/malva/examples/malva_tools/3M-february-2018.txt'
    state.update(source=url,expected_HTTP_content_length=115512960,reference_provenance='Public MALVA distribution copy; compare fingerprint to simpleaf registry when available')
    for attempt in range(1,6):
        if p.exists() and p.stat().st_size==115512960:break
        with (D/'reference515.curl.log').open('ab') as log:
            r=subprocess.run(['curl','-q','--proxy','','--noproxy','*','--fail','--location','--continue-at','-','--connect-timeout','30','--max-time','7200','--speed-time','180','--speed-limit','1024','--output',str(p),url],stdout=log,stderr=log)
        state.update(attempt=attempt,curl_exit=r.returncode,bytes=p.stat().st_size if p.exists() else 0);save()
        if r.returncode==0:break
        time.sleep(30)
    assert p.stat().st_size==115512960,'Reference incomplete; do not evaluate partial reference'
    state['reference_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    registry=json.loads((D/'simpleaf_chemistries.json').read_text(encoding='utf-8'))['10xv3']
    state['simpleaf_registry_fingerprint']=registry['plist_name']
    state['SHA256_equals_registry_fingerprint']=state['reference_sha256']==registry['plist_name']
    state['stage']='barcode_comparison';save()
    whitelist=set()
    with p.open(encoding='ascii') as f:
        for line in f:
            seq=line.strip();assert len(seq)==16 and set(seq)<=set('ACGT');whitelist.add(seq)
    state['reference_unique_barcodes']=len(whitelist)
    comp=str.maketrans('ACGT','TGCA');counts=collections.Counter();top=[]
    with (Q/'alevin/raw_cb_frequency.txt').open(encoding='ascii') as f:
        for line in f:
            cb,number=line.split();number=int(number)
            match=cb in whitelist;rc=cb.translate(comp)[::-1] in whitelist
            counts['observed_barcodes']+=1;counts['observed_reads']+=number
            counts['exact_barcodes']+=match;counts['exact_reads']+=number*match
            counts['reverse_complement_barcodes']+=rc;counts['reverse_complement_reads']+=number*rc
            if len(top)<20:top.append(dict(barcode=cb,reads=number,exact_reference_match=match,RC_reference_match=rc))
    called=(Q/'alevin/whitelist.txt').read_text(encoding='ascii').splitlines()
    with gzip.open(B/'results/pilot_cell_qualification_509/cell_qualification.tsv.gz','rt',encoding='utf-8') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    good=[r for r in rows if r['barcode'] in set(called) and r['QC_pass']=='True']
    state.update(counts=dict(counts),top_barcodes=top,final_called=len(called),final_called_exact_reference_match=sum(cb in whitelist for cb in called),QC_pass_called=len(good),QC_pass_called_exact_reference_match=sum(r['barcode'] in whitelist for r in good),stage='complete',finished=now())
except Exception:state.update(stage='failed',error=traceback.format_exc(),finished=now())
save();print(json.dumps(state),flush=True)

