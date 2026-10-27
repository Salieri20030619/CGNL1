from pathlib import Path
import os,select,json,gzip,hashlib,subprocess,datetime,shutil
O=Path.cwd()/'mechanism_investigation_v2';R=O/'reference/gencode_v47';receipt=O/'metadata/REFERENCE_DOWNLOAD_COMPLETE.json';install=O/'metadata/QUANTIFIER_INSTALL.json'
def save(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
# Native process-exit wait; not periodic agent-level polling.
if not receipt.exists():
 pids=[]
 for d in Path('/proc').iterdir():
  if d.name.isdigit():
   try:
    args=(d/'cmdline').read_bytes().split(b'\0')
    if any(a.endswith(b'/cgnl_mech41_acquire_reference.py') for a in args):pids.append(int(d.name))
   except (FileNotFoundError,PermissionError):pass
 assert len(pids)==1,('reference downloader not uniquely live',pids)
 fd=os.pidfd_open(pids[0]);save(O/'metadata/INDEX_DEPENDENCY_WAIT.json',{'pid':os.getpid(),'waiting_on_reference_pid':pids[0],'method':'Linux pidfd process exit event','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});select.select([fd],[],[]);os.close(fd)
assert receipt.exists(),'Downloader exited without completion receipt';rs=json.loads(receipt.read_text());assert all(r.get('success') for r in rs),rs
q=json.loads(install.read_text());assert q['exit_code']==0 and q['version']=='salmon 1.10.3'
md={}
for line in (R/'MD5SUMS').read_text(encoding='utf-8').splitlines():
 parts=line.split()
 if len(parts)==2:md[parts[1].lstrip('*').removeprefix('./')]=parts[0]
verified=[]
for n in ['gencode.v47.transcripts.fa.gz','GRCh38.primary_assembly.genome.fa.gz']:
 h=hashlib.md5()
 with (R/n).open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 assert n in md and h.hexdigest()==md[n],n;verified.append({'name':n,'md5':h.hexdigest()})
gentrome=R/'gentrome.fa';decoys=R/'decoys.txt';mapping=R/'transcript_to_gene.tsv';assert not gentrome.exists()
with gentrome.open('wb') as out,mapping.open('w',encoding='utf-8',newline='') as mp:
 mp.write('transcript_id\tgene_id\tgene_symbol\n')
 with gzip.open(R/'gencode.v47.transcripts.fa.gz','rb') as f:
  for line in f:
   if line.startswith(b'>'):
    parts=line[1:].decode('utf-8').strip().split('|');assert len(parts)>=6;mp.write(parts[0]+'\t'+parts[1]+'\t'+parts[5]+'\n')
   out.write(line)
 with gzip.open(R/'GRCh38.primary_assembly.genome.fa.gz','rb') as f,decoys.open('w',encoding='utf-8') as dc:
  for line in f:
   if line.startswith(b'>'):dc.write(line[1:].split()[0].decode('utf-8')+'\n')
   out.write(line)
exe=O/'software/salmon_1_10_3/bin/salmon';idx=R/'salmon_index';cmd=[str(exe),'index','-t',str(gentrome),'-d',str(decoys),'-i',str(idx),'--gencode','-k','31','-p','8']
save(O/'metadata/INDEX_SPECIFICATION.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs':verified,'command':cmd,'pid':os.getpid(),'reference':'GENCODE47 allchromosomaltranscripts withGRCh38primarygenomedecoy','threads':8})
with (O/'logs/salmon_index.log').open('wb') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
assert r.returncode==0,('index exit',r.returncode)
assert (idx/'versionInfo.json').exists() and (idx/'info.json').exists()
save(O/'metadata/INDEX_COMPLETE.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':r.returncode,'index':str(idx),'inputs':verified,'versionInfo':json.loads((idx/'versionInfo.json').read_text()),'info':json.loads((idx/'info.json').read_text())});print('Index complete',str(idx),flush=True)
