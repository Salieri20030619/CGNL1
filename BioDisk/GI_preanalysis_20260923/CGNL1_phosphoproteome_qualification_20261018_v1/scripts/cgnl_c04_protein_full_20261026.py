from pathlib import Path
import json,hashlib,subprocess,datetime,os,fcntl,xml.etree.ElementTree as E,collections,traceback,shutil,csv,gzip,math
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_totalprotein_full_20261026';P=O/'C04_totalprotein_pilot_20261026';S=O/'C04_search_pilot'
lock=(W/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for b in iter(lambda:stream.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(path,x):
 tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(path)
state=dict(pid=os.getpid(),started=now(),files={},complete=False,whole_project_complete=False)
def update(name,stage,**kwargs):
 state['files'][name]=dict(stage=stage,utc=now(),**kwargs);state['current_file']=name;save(W/'STATUS.json',state);print(json.dumps(dict(file=name,stage=stage,**kwargs)),flush=True)
def validate_outputs(folder):
 r=json.loads((folder/'SEARCH_EXIT.json').read_text());assert r['exit']==0
 assert {'psms.txt','psms.pin','psms.pep.xml'}<={x['name'] for x in r['files']}
 for x in r['files']:
  p=folder/x['name'];assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256']
 return r
try:
 assert not (W/'EXIT.json').exists(), 'Prior terminal state requires audit/versioned recovery'
 spec=json.loads((W/'FULL_SCOPE_SPEC_BEFORE_NEW_RESULTS.json').read_text());assert len(spec['files'])==24
 gate=P/'TWO_BATCH_PILOT_ACCEPTANCE.json';assert sha(gate)==spec['pilot_acceptance_sha256'];assert json.loads(gate.read_text())['two_batch_acquisition_and_search_pilot_pass']
 template=Path(spec['parameter_template']);assert sha(template)==spec['parameter_sha256'];assert sha(S/'mouse_plus_contaminants.fasta')==spec['database_sha256']
 search_spec=json.loads((P/'protein_set1_fraction_1/search_CID/SPEC_BEFORE_SEARCH.json').read_text());assert sha(S/'comet.linux.exe')==search_spec['engine_sha256']
 acquisition_spec=json.loads((P/'SPEC_BEFORE_RESULTS.json').read_text());assert sha(O/'software/ThermoRawFileParser_2_dev/ThermoRawFileParser')==acquisition_spec['converter_sha256']
 assert shutil.disk_usage(W).free>40*1024**3
 for item in spec['files']:
  name=item['name'];stem=Path(name).stem;folder=W/stem;accept=W/(name+'.ACCEPTANCE.json')
  if name in spec['reuse_verified_pilot_files']:
   src=P/stem;a=json.loads((P/(name+'.ACCEPTANCE.json')).read_text());assert a['acquisition_pass'];assert sha(src/name)==a['raw_sha256'];assert sha(src/(stem+'.mzML'))==a['mzML_sha256'];assert sha(src/'search_CID/comet.params')==spec['parameter_sha256'];validate_outputs(src/'search_CID')
   assert not folder.exists();folder.symlink_to(src,target_is_directory=True);save(accept,a);update(name,'reused_verified_pilot',source=str(src));continue
  folder.mkdir(exist_ok=True);raw=folder/name;part=folder/(name+'.partial');assert not raw.exists(),'Untracked existing RAW requires audit'
  assert shutil.disk_usage(W).free>15*1024**3;update(name,'downloading',expected_bytes=item['bytes'],prefix_bytes=part.stat().st_size if part.exists() else 0)
  cmd=['curl','-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time','43200','--speed-limit','1024','--speed-time','180','--silent','--show-error']
  if part.exists():cmd+=['--continue-at','-']
  cmd+=['--output',str(part),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name]
  r=subprocess.run(cmd,capture_output=True,text=True);save(folder/'TRANSFER.json',dict(exit=r.returncode,bytes=part.stat().st_size if part.exists() else 0,stderr=r.stderr,utc=now(),TLS_verification_enabled=True))
  assert r.returncode==0 and part.stat().st_size==item['bytes'], 'Transfer incomplete; preserve prefix and audit before resuming'
  with part.open('rb') as stream:assert stream.read(2)==b'\x01\xa1'
  digest=sha(part);part.rename(raw);save(folder/'RAW_RECEIPT.json',dict(bytes=item['bytes'],sha256=digest,utc=now()));update(name,'converting',raw_sha256=digest)
  converter=O/'software/ThermoRawFileParser_2_dev/ThermoRawFileParser'
  with (folder/'conversion.log').open('w',encoding='utf-8') as log:r=subprocess.run([str(converter),'-i='+str(raw),'-o='+str(folder),'-f=2','-m=0'],stdout=log,stderr=subprocess.STDOUT)
  save(folder/'CONVERSION_EXIT.json',dict(exit=r.returncode,utc=now()));assert r.returncode==0
  mz=folder/(stem+'.mzML');ns={'m':'http://psi.hupo.org/ms/mzml'};counts=collections.Counter();modes=collections.Counter();ids=set();missing=0;nr=0
  for _,node in E.iterparse(mz,events=['end']):
   if node.tag!='{'+ns['m']+'}spectrum':continue
   cv={x.attrib['accession']:x.attrib.get('value') for x in node.findall('m:cvParam',ns)};level=int(cv['MS:1000511']);counts[level]+=1
   if level==2:
    for x in node.findall('.//m:activation/m:cvParam',ns):
     if 'dissociation' in x.attrib.get('name',''):modes[x.attrib['name']]+=1
   if level==3:
    refs=[x.attrib.get('spectrumRef') for x in node.findall('m:precursorList/m:precursor',ns)];nr+=len(refs);missing+=sum(x not in ids for x in refs)
   ids.add(node.attrib['id']);node.clear()
  assert all(counts[x]>0 for x in [1,2,3]) and missing==0 and nr>0
  assert dict(modes)=={'collision-induced dissociation':counts[2]},'Unexpected acquisition mode; do not apply wrong search'
  a=dict(name=name,utc=now(),raw_sha256=digest,mzML_bytes=mz.stat().st_size,mzML_sha256=sha(mz),counts=dict(counts),MS2_activation_counts=dict(modes),MS3_references=nr,unresolved_MS3_references=missing,conversion_exit=0,acquisition_pass=True);save(accept,a)
  out=folder/'search_CID';out.mkdir(exist_ok=True);assert not (out/'RUN.json').exists();assert sha(template)==spec['parameter_sha256'];(out/'comet.params').write_bytes(template.read_bytes());(out/'input.mzML').symlink_to(mz)
  available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')));assert available>16*1024**2
  save(out/'INPUT_SPEC.json',dict(file=name,params_sha256=spec['parameter_sha256'],mzML_sha256=a['mzML_sha256'],database_sha256=spec['database_sha256'],utc=now()));update(name,'searching_CID')
  with (out/'search.log').open('w',encoding='utf-8') as log:
   child=subprocess.Popen([str(S/'comet.linux.exe'),'-P'+str(out/'comet.params'),'-N'+str(out/'psms'),str(out/'input.mzML')],cwd=out,stdout=log,stderr=subprocess.STDOUT);save(out/'RUN.json',dict(pid=child.pid,parent_pid=os.getpid(),utc=now()));code=child.wait()
  files=[dict(name=x.name,bytes=x.stat().st_size,sha256=sha(x)) for x in out.glob('psms*') if x.is_file()];save(out/'SEARCH_EXIT.json',dict(exit=code,files=files,utc=now()));validate_outputs(out);update(name,'acquisition_and_search_complete')
 assert len(state['files'])==24
 A=W/'joint_identification';A.mkdir(exist_ok=True);inputs=[];rows=[];seen=set();csv.field_size_limit(20000000)
 for item in spec['files']:
  out=W/Path(item['name']).stem/'search_CID';validate_outputs(out);assert sha(out/'comet.params')==spec['parameter_sha256'];f=out/'psms.txt';inputs.append(dict(file=item['name'],path=str(f),sha256=sha(f)))
  with f.open(encoding='utf-8') as stream:
   assert stream.readline().startswith('CometVersion 2026.03')
   for x in csv.DictReader(stream,delimiter='\t'):
    assert int(x['num'])==1;key=(item['name'],x['scan'],x['charge']);assert key not in seen;seen.add(key);assert '79.966' not in x['modifications'];e=float(x['e-value']);assert math.isfinite(e) and e>=0
    rows.append(dict(file=item['name'],scan=int(x['scan']),charge=int(x['charge']),e_value=e,is_decoy='DECOY_' in x['protein'],is_contaminant='CON__' in x['protein'],peptide=x['plain_peptide'],protein=x['protein']))
 save(A/'INPUT_MANIFEST.json',dict(utc=now(),inputs=inputs,complete24scope=True,parameter_sha256=spec['parameter_sha256']))
 groups=collections.defaultdict(list)
 for x in rows:groups[x['e_value']].append(x)
 targets=decoys=0;ordered=[]
 for e,g in sorted(groups.items()):
  decoys+=sum(x['is_decoy'] for x in g);targets+=sum(not x['is_decoy'] for x in g);ordered.append((g,min(1,(decoys+1)/max(targets,1))))
 q=1
 for g,fdr in reversed(ordered):
  q=min(q,fdr)
  for x in g:x['q_fullscope']=q
 with gzip.open(A/'PSM_QVALUES.tsv.gz','wt',encoding='utf-8',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
 geneids=set()
 for line in (S/'mouse_plus_contaminants.fasta').open():
  if line.startswith('>') and 'GN=Cgnl1 ' in line:geneids.add(line[1:].split()[0])
 cg=[x for x in rows if not x['is_decoy'] and geneids.intersection(x['protein'].split(','))];save(A/'CGNL1_CANDIDATES.json',dict(records=cg))
 save(A/'RESULTS.json',dict(utc=now(),PSMs=len(rows),targets=targets,decoys=decoys,passing_noncontaminant_PSMs=sum(not x['is_decoy'] and not x['is_contaminant'] and x['q_fullscope']<=.01 for x in rows),CGNL1_candidate_PSMs=len(cg),independent_QA_pending=True,protein_FDR_not_established=True,protein_quantification_pending=True,whole_project_complete=False))
 state.update(complete=True,stage='full24_acquisition_search_jointPSM_complete_QA_pending')
except Exception:
 state.update(stage='failed',error=traceback.format_exc());raise
finally:
 state['finished']=now();save(W/'STATUS.json',state);save(W/'EXIT.json',dict(complete=state['complete'],error=state.get('error'),utc=now(),whole_project_complete=False));print(json.dumps(state),flush=True)
