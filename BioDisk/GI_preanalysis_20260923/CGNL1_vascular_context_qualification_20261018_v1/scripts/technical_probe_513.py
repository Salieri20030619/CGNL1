from pathlib import Path
import collections,datetime,gzip,hashlib,json,os,subprocess,traceback
import numpy as np
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
root=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018');D=root/'technical_probe_513'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(data):
    p=B/'metadata/TECHNICAL_PROBE_513.json';p.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
state=dict(pid=os.getpid(),started=now(),stage='sampling',scope='first 1000000 read pairs; technical layout probe only, not whole-library quality estimate',biological_inference=False)
save(state)
try:
    whitelist_path=D/'3M-february-2018.txt.gz'
    whitelist=None
    if whitelist_path.exists():
        with gzip.open(whitelist_path,'rt',encoding='ascii') as f:
            whitelist={line.split()[0] for line in f if line.strip()}
        assert len(whitelist)>1000000 and all(len(x)==16 and set(x)<=set('ACGT') for x in whitelist)
        state['whitelist']=dict(sha256=hashlib.sha256(whitelist_path.read_bytes()).hexdigest(),barcodes=len(whitelist))
    else:
        state['whitelist']=dict(available=False,matching_performed=False,reason='Original vendor URL absent; developer-maintained simpleaf registry located, public Box download timed out. Continue independent quality and strand probes.')
    a=root/'r1_fresh_495/PT55_TU_S1_L002_R1_001.fastq.gz';c=root/'proxy_redownload_492/PT55_TU_S1_L002_R2_001.fastq.gz'
    output=D/'R2_prefix_1M.fastq.gz';assert not output.exists()
    quality=[np.zeros(28),np.zeros(91)];q30=[np.zeros(28),np.zeros(91)];base_counts=[np.zeros((28,5),dtype=np.int64),np.zeros((91,5),dtype=np.int64)]
    hits=collections.Counter();bc=collections.Counter();comp=str.maketrans('ACGTN','TGCAN');n=0
    with gzip.open(a,'rt',encoding='ascii') as f1,gzip.open(c,'rt',encoding='ascii') as f2,gzip.open(output,'wt',encoding='ascii',compresslevel=1) as out:
        for k in range(1000000):
            r1=[f1.readline() for _ in range(4)];r2=[f2.readline() for _ in range(4)]
            if not r1[0]: break
            assert all(r1) and all(r2) and r1[0].split()[0]==r2[0].split()[0]
            s1=r1[1].strip();s2=r2[1].strip();assert len(s1)==28 and len(s2)==91
            out.writelines(r2);bc[s1[:16]]+=1
            candidates={'first16':s1[:16],'last16':s1[-16:],'RC_first16':s1.translate(comp)[::-1][:16],'RC_of_first16':s1[:16].translate(comp)[::-1],'R2_first16':s2[:16]}
            if whitelist is not None:
                for name,seq in candidates.items(): hits[name]+=seq in whitelist
            for i,(seq,record) in enumerate([(s1,r1),(s2,r2)]):
                qq=np.frombuffer(record[3].strip().encode('ascii'),dtype=np.uint8).astype(np.int16)-33
                assert len(qq)==len(seq) and (qq>=0).all()
                quality[i]+=qq;q30[i]+=qq>=30
                for code,base in enumerate('ACGTN'):base_counts[i][:,code]+=np.frombuffer(seq.encode('ascii'),dtype=np.uint8)==ord(base)
            n+=1
    state.update(sampled_pairs=n,barcode_exact_match_fractions={key:value/n for key,value in hits.items()} if whitelist is not None else None,top_sample_barcodes=[dict(barcode=x,reads=count,in_whitelist=(x in whitelist) if whitelist is not None else None) for x,count in bc.most_common(10)],mean_quality_by_position=[(x/n).tolist() for x in quality],Q30_fraction_by_position=[(x/n).tolist() for x in q30],base_counts=[x.tolist() for x in base_counts],stage='auto_strand_mapping')
    save(state)
    spec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_496.json').read_text(encoding='utf-8'))
    exe=Path(spec['reference']).parents[2]/'software/salmon_1_10_3/bin/salmon'
    cmd=[str(exe),'quant','-i',spec['reference'],'-l','A','-r',str(output),'-p','8','-o',str(D/'auto_strand_quant'),'--validateMappings']
    state['command']=cmd;save(state)
    with (D/'auto_strand.log').open('wb') as log:r=subprocess.run(cmd,stdout=log,stderr=log)
    state['mapping_exit']=r.returncode;assert r.returncode==0
    state['mapping_meta']=json.loads((D/'auto_strand_quant/aux_info/meta_info.json').read_text(encoding='utf-8'))
    state['library_formats']=json.loads((D/'auto_strand_quant/lib_format_counts.json').read_text(encoding='utf-8'))
    state['stage']='complete'
except Exception:state.update(stage='failed',error=traceback.format_exc())
state['finished']=now();save(state)
print(json.dumps({k:v for k,v in state.items() if k not in ['mean_quality_by_position','Q30_fraction_by_position','base_counts']}),flush=True)

