from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,datetime,subprocess,hashlib
import pandas as pd,numpy as np
O=Path.cwd()/'mechanism_investigation_v2';plan=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text());runs={}
for t in plan['tasks']:runs.setdefault(t['run'],[]).append(t)
assert len(runs)==12 and all(len(v)==2 for v in runs.values());Q=O/'processed/DRA014736_quant';Q.mkdir(parents=True,exist_ok=True);done={};fail={}
def save(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def live(script):
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    if any(a.endswith(('/'+script).encode()) for a in (p/'cmdline').read_bytes().split(b'\0')):return True
   except (FileNotFoundError,PermissionError):pass
 return False
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'salmon':'1.10.3','reference':'GENCODE47 transcriptome +GRCh38primarydecoys,authorMD5required','library_type':'automatic, inspect inferred strand for all samples','options':['--validateMappings','--seqBias','--gcBias'],'threads':8,'concurrency':1,'inputs':'only pairedFASTQwithENAmd5verified','output':'fulltranscript quant and gene sumTPM/estimatedNumReads, notUMI','statistical_role':'3donors descriptive response, no inferential claim from readcounts or number oflibraries','monitor':'single persistent900secondwait when prerequisites absent'}
save(O/'DRA_QUANTIFICATION_SPECIFICATION_v2.json',spec)
while len(done)+len(fail)<len(runs):
 progressed=False;index=O/'metadata/INDEX_COMPLETE.json';receipts=O/'metadata/DRA_DOWNLOAD_RECEIPTS.json'
 if index.exists() and receipts.exists():
  idx=json.loads(index.read_text());assert idx['exit_code']==0;rr=json.loads(receipts.read_text());verified={r['name']:r for r in rr if r.get('success')};mapping=pd.read_csv(O/'reference/gencode_v47/transcript_to_gene.tsv',sep='\t',dtype=str);assert mapping.transcript_id.is_unique
  for run,files in runs.items():
   if run in done or run in fail or not all(t['name'] in verified for t in files):continue
   files=sorted(files,key=lambda t:t['name']);paths=[O/'raw/DRA014736'/t['name'] for t in files];dest=Q/run
   assert not dest.exists(),'existing result requires review'
   cmd=[str(O/'software/salmon_1_10_3/bin/salmon'),'quant','-i',idx['index'],'-l','A','-1',str(paths[0]),'-2',str(paths[1]),'--validateMappings','--seqBias','--gcBias','-p','8','-o',str(dest)]
   item={'run':run,'donor':files[0]['donor'],'condition':files[0]['condition'],'command':cmd,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
   for t,p in zip(files,paths):
    h=hashlib.md5()
    with p.open('rb') as f:
     while b:=f.read(4*1024*1024):h.update(b)
    assert p.stat().st_size==t['bytes'] and h.hexdigest()==t['md5']
   with (O/'logs'/('salmon_'+run+'.log')).open('wb') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   item['exit_code']=r.returncode
   if r.returncode:fail[run]=item
   else:
    a=pd.read_csv(dest/'quant.sf',sep='\t');assert np.isfinite(a[['TPM','NumReads']].values).all() and (a[['TPM','NumReads']]>=0).all().all();v=a.merge(mapping,left_on='Name',right_on='transcript_id',how='left',validate='one_to_one');assert v.gene_id.notna().all();g=v.groupby(['gene_id','gene_symbol'])[['TPM','NumReads']].sum();g.to_csv(dest/'gene_quant.tsv.gz',sep='\t',encoding='utf-8');meta=json.loads((dest/'aux_info/meta_info.json').read_text());item.update(meta_info=meta,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());done[run]=item
   save(O/'metadata/DRA_QUANTIFICATION_STATUS.json',{'pid':os.getpid(),'done':done,'fail':fail,'expected':12});print('QUANT',run,'exit',r.returncode,flush=True);progressed=True
 if len(done)+len(fail)==12:break
 if not live('cgnl_mech67_resumable_transfers.py') and not live('cgnl_mech45c_index_when_ready.py'):
  if not progressed:break
 if not progressed:time.sleep(900)
save(O/'metadata/DRA_QUANTIFICATION_COMPLETE.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':len(done)==12 and not fail,'done_runs':list(done),'failed_runs':fail,'missing_runs':[r for r in runs if r not in done and r not in fail]});print('Quantification controller finished',len(done),'of12',flush=True)
