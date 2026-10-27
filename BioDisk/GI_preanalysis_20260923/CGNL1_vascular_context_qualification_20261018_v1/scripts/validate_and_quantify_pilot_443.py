from pathlib import Path
import os,select,json,gzip,collections,time,subprocess,datetime
b=Path(__file__).resolve().parents[1];spec=json.loads((b/'metadata/PILOT_QUANTIFICATION_SPEC_443.json').read_text());launch=json.loads(Path(spec['download_dependency']).read_text());receipt=Path(launch['receipt'])
if not receipt.exists():
 try:
  fd=os.pidfd_open(launch['wrapper_pid']);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
 except ProcessLookupError:pass
assert receipt.exists(), 'Download process ended without exit receipt'
assert json.loads(receipt.read_text())['exit_code']==0, 'Download failed; do not quantify'
d=json.loads((b/'metadata/PILOT_DOWNLOAD_COMPLETE_441.json').read_text());assert d['all_md5_verified']
files={('R1' if '_R1_' in f['name'] else 'R2'):Path(f['local_path']) for f in d['files']};manifest=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text())
for item in manifest['items']:
 p=Path(manifest['destination'])/item['name'];assert p.stat().st_size==item['bytes']
counts=[collections.Counter(),collections.Counter()];n=0;start=time.monotonic()
with gzip.open(files['R1'],'rt',encoding='ascii',newline='') as a,gzip.open(files['R2'],'rt',encoding='ascii',newline='') as c:
 while True:
  h1=a.readline();h2=c.readline()
  if not h1 or not h2:
   assert not h1 and not h2,'unequal FASTQ lengths';break
  s1,p1,q1=a.readline().rstrip('\r\n'),a.readline(),a.readline().rstrip('\r\n')
  s2,p2,q2=c.readline().rstrip('\r\n'),c.readline(),c.readline().rstrip('\r\n')
  assert h1.startswith('@') and h2.startswith('@') and p1.startswith('+') and p2.startswith('+')
  assert len(s1)==len(q1) and len(s2)==len(q2) and len(s1)==28 and len(s2)==91
  assert h1.split()[0].removesuffix('/1')==h2.split()[0].removesuffix('/2')
  counts[0][len(s1)]+=1;counts[1][len(s2)]+=1;n+=1
  if n%1000000==0:
   (b/'metadata/FASTQ_QA_PROGRESS_443.json').write_text(json.dumps({'reads':n,'elapsed_seconds':time.monotonic()-start,'pid':os.getpid()}),encoding='utf-8')
assert n>0
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':True,'read_pairs':n,'R1_lengths':dict(counts[0]),'R2_lengths':dict(counts[1]),'gzip_full_EOF_CRC_checked':True,'paired_identifiers_all_matched':True}
(b/'metadata/FASTQ_QA_COMPLETE_443.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
exe=Path(spec['reference']).parents[2]/'software/salmon_1_10_3/bin/salmon'
cmd=[str(exe),'alevin','-l','ISR','-1',str(files['R1']),'-2',str(files['R2']),'--chromiumV3','-i',spec['reference'],'-p',str(spec['threads']),'-o',spec['output'],'--tgMap',spec['transcript_gene_map'],'--dumpMtx','--dumpFeatures']
(b/'metadata/ALEVIN_COMMAND_443.json').write_text(json.dumps(cmd,indent=2),encoding='utf-8')
print('QUANT_START',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
r=subprocess.run(cmd);assert r.returncode==0, 'Alevin nonzero exit'
for name in spec['required_outputs']:
 p=Path(spec['output'])/name;assert p.is_file() and p.stat().st_size>0,name
(b/'metadata/PILOT_QUANT_COMPLETE_443.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':0,'required_files_present':True,'biological_validation':False,'output':spec['output']},indent=2),encoding='utf-8')
