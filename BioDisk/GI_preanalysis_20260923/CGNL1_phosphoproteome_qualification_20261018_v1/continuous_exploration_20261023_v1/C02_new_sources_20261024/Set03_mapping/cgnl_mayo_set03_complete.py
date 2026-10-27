from pathlib import Path
import os,select,json,subprocess,hashlib,csv,math,statistics,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024');O=P/'Set03_mapping'
pid=474238;target=str(P/'Set03_Phospho/evidence.txt.deflate.partial')
if Path(f'/proc/{pid}').exists():
    fds=Path(f'/proc/{pid}/fd');links=[]
    for f in fds.iterdir():
        try:links.append(str(f.readlink()))
        except FileNotFoundError:pass
    assert target in links,'PID no longer owns expected partial; inspect before adoption'
    fd=os.pidfd_open(pid);poll=select.poll();poll.register(fd,select.POLLIN)
    (O/'COMPLETION_WAITER.json').write_text(json.dumps({'pid':os.getpid(),'adopted_download_pid':pid,'native_pidfd_wait':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2),encoding='utf-8')
    print('Waiting on verified download PID via pidfd; no periodic polling',flush=True)
    poll.poll();os.close(fd)
spec=O/'DESCRIPTIVE_SPEC_BEFORE_RESULTS.json'
assert hashlib.sha256(spec.read_bytes()).hexdigest()=='e3773c094398f28f8e284a4c80b2d29595c2474c95e181452031ba572afd5079'
run=subprocess.run(['python3',str(O/'evidence_gate.py')],capture_output=True,text=True)
(O/'EVIDENCE_GATE_RUN.json').write_text(json.dumps({'exit':run.returncode,'stdout':run.stdout,'stderr':run.stderr},indent=2),encoding='utf-8')
assert run.returncode==0,run.stderr or run.stdout
g=json.loads((P/'Set03_EVIDENCE_GATE.json').read_text())
if not g['pass']:
    (O/'DESCRIPTIVE_NON_EVALUABLE.json').write_text(json.dumps({'reason':'frozen evidence gate failed','gate_path':str(P/'Set03_EVIDENCE_GATE.json'),'no_biological_tests':True},indent=2),encoding='utf-8')
    print('Evidence gate failed; no contrast computed',flush=True);raise SystemExit(0)
assert json.loads((O/'CHANNEL_INDEPENDENT_QA.json').read_text())['pass']
assert json.loads((O/'NORMAL_MAPPING_INDEPENDENT_QA.json').read_text())['pass']
output=O/'DESCRIPTIVE_RESULTS.json';assert not output.exists(),'Preserve prior result'
rows=[]
for sample,channel,group in [('iCCA27',10,'tumor'),('iCCA29',12,'tumor'),('iCCA30',13,'tumor'),('Normal5',14,'normal'),('Normal6',15,'normal')]:
    a=g['sites']['486']['site_row_mono_reporters'][channel-1];b=g['sites']['283']['site_row_mono_reporters'][channel-1]
    assert a>0 and b>0 and math.isfinite(a) and math.isfinite(b)
    rows.append({'sample':sample,'channel':channel,'group':group,'S486':a,'S283':b,'log2ratio':math.log2(a/b)})
with (O/'DESCRIPTIVE_INPUT.csv').open('w',encoding='utf-8',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
comparisons=[{'tumor':a['sample'],'normal':b['sample'],'delta':a['log2ratio']-b['log2ratio']} for a in rows[:3] for b in rows[3:]]
rs="""args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(file.path(p,'DESCRIPTIVE_INPUT.csv'));r<-log2(d$S486)-log2(d$S283)
stopifnot(max(abs(r-d$log2ratio))<1e-12)
out<-do.call(rbind,lapply(1:3,function(i)data.frame(tumor=d$sample[i],normal=d$sample[4:5],delta=r[i]-r[4:5])))
write.csv(out,file.path(p,'DESCRIPTIVE_R.csv'),row.names=FALSE)
"""
(O/'DESCRIPTIVE_INDEPENDENT.R').write_text(rs,encoding='utf-8')
run=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'DESCRIPTIVE_INDEPENDENT.R'),str(O)],capture_output=True,text=True);assert run.returncode==0,run.stderr
with (O/'DESCRIPTIVE_R.csv').open() as h:rr=list(csv.DictReader(h))
assert len(rr)==6;errors=[]
for a,b in zip(comparisons,rr):
    assert a['tumor']==b['tumor'] and a['normal']==b['normal'];errors.append(abs(a['delta']-float(b['delta'])))
assert max(errors)<1e-12
per=[]
for row in rows[:3]:
    dd=[x['delta'] for x in comparisons if x['tumor']==row['sample']];per.append({'tumor':row['sample'],'median':statistics.median(dd),'range':[min(dd),max(dd)]})
result={'spec_sha256':hashlib.sha256(spec.read_bytes()).hexdigest(),'samples':rows,'comparisons':comparisons,'per_tumor':per,'tumor_n':3,'normal_specimens':2,'no_pairing_inferred':True,'P':None,'FDR':None,'R_QA_pass':True,'R_max_error':max(errors),'no_new_biological_tests':True,'whole_project_complete':False}
output.write_text(json.dumps(result,indent=2),encoding='utf-8');assert json.loads(output.read_text())==result
print(json.dumps({'per_tumor':per,'R_QA_pass':True}),flush=True)
