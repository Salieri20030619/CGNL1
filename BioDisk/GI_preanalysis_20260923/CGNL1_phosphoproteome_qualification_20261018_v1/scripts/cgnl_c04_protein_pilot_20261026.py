from pathlib import Path
import json,hashlib,subprocess,datetime,os,fcntl,xml.etree.ElementTree as E,collections,traceback,shutil
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1')
W=P/'C04_totalprotein_pilot_20261026';W.mkdir(exist_ok=True)
lock=(W/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def save(n,x):
 x.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),whole_project_complete=False)
 tmp=W/(n+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8');tmp.replace(W/n)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(8388608),b''):h.update(b)
 return h.hexdigest()
state=dict(pid=os.getpid(),files={},stage='starting',complete=False)
try:
 assert not (W/'EXIT.json').exists(), 'Previous terminal receipt requires audit'
 assert shutil.disk_usage(W).free>15*1024**3
 source=P/'C04_sources/PAK_file_manifest.json';allrows=json.loads(source.read_text())['row_data']
 items=[dict(name=n,bytes=next(x['size'] for x in allrows if x['name']==n)) for n in ['protein_set1_fraction_1.raw','protein_set2_fraction_1.raw']]
 binary=P/'software/ThermoRawFileParser_2_dev/ThermoRawFileParser'
 save('SPEC_BEFORE_RESULTS.json',dict(files=items,selection='fraction1 from each protein TMT batch, fixed before protein results',source_manifest_sha256=sha(source),converter_sha256=sha(binary),read_gate=['RAW exact source bytes and magic','conversion exit0','MS1 MS2 MS3 present','all MS3 precursor references resolve'],search_started=False,full24_expansion_requires_pilot_acceptance=True,CGNL1_absence_not_inferred_from_pilot=True,protein_phospho_channels_distinct=True,disease_tests=0))
 for item in items:
  name=item['name'];d=W/Path(name).stem;d.mkdir(exist_ok=True);raw=d/name;part=d/(name+'.partial');assert not raw.exists()
  state.update(stage='downloading',current_file=name);state['files'][name]={'stage':'downloading'};save('STATUS.json',state)
  cmd=['curl','-q','--ipv4','--noproxy','*','--ssl-reqd','--tls-max','1.2','--connect-timeout','20','--max-time','43200','--speed-limit','1024','--speed-time','180','--silent','--show-error']
  if part.exists():cmd+=['--continue-at','-']
  cmd+=['--output',str(part),'ftp://massive-ftp.ucsd.edu/v01/MSV000081198/raw/'+name]
  r=subprocess.run(cmd,capture_output=True,text=True);save(name+'.TRANSFER.json',dict(exit=r.returncode,bytes=part.stat().st_size if part.exists() else 0,stderr=r.stderr,TLS_verification_enabled=True));assert r.returncode==0 and part.stat().st_size==item['bytes']
  with part.open('rb') as fh:assert fh.read(2)==b'\x01\xa1'
  digest=sha(part);part.rename(raw);save(name+'.RAW_RECEIPT.json',dict(bytes=raw.stat().st_size,sha256=digest))
  state['files'][name]={'stage':'converting','raw_sha256':digest};state['stage']='converting';save('STATUS.json',state)
  with (d/'conversion.log').open('w',encoding='utf-8') as log:r=subprocess.run([str(binary),'-i='+str(raw),'-o='+str(d),'-f=2','-m=0'],stdout=log,stderr=subprocess.STDOUT)
  save(name+'.CONVERSION_EXIT.json',dict(exit=r.returncode));assert r.returncode==0
  mzml=d/(Path(name).stem+'.mzML');assert mzml.exists()
  ns={'m':'http://psi.hupo.org/ms/mzml'};counts=collections.Counter();modes=collections.Counter();ids=set();missing=[];refs_count=0
  for _,node in E.iterparse(mzml,events=['end']):
   if node.tag!='{'+ns['m']+'}spectrum':continue
   cv={x.attrib['accession']:x.attrib.get('value') for x in node.findall('m:cvParam',ns)};level=int(cv['MS:1000511']);counts[level]+=1
   if level==2:
    for x in node.findall('.//m:activation/m:cvParam',ns):
     if 'dissociation' in x.attrib.get('name',''):modes[x.attrib['name']]+=1
   if level==3:
    refs=[x.attrib.get('spectrumRef') for x in node.findall('m:precursorList/m:precursor',ns)];refs_count+=len(refs);missing.extend(x for x in refs if x not in ids)
   ids.add(node.attrib['id']);node.clear()
  assert all(counts[k]>0 for k in [1,2,3]) and not missing and refs_count>0
  accept=dict(name=name,raw_sha256=digest,mzML_bytes=mzml.stat().st_size,mzML_sha256=sha(mzml),counts=dict(counts),MS2_activation_counts=dict(modes),MS3_references=refs_count,unresolved_MS3_references=len(missing),conversion_exit=0,acquisition_pass=True,protein_identification_not_run=True)
  save(name+'.ACCEPTANCE.json',accept);state['files'][name]={'stage':'acquisition_pass','receipt':str(W/(name+'.ACCEPTANCE.json'))};save('STATUS.json',state);print(json.dumps(accept),flush=True)
 state.update(stage='pilot_acquisition_complete_search_pending',complete=True);save('STATUS.json',state);save('EXIT.json',dict(acquisition_complete=True,protein_search_complete=False))
except Exception:
 state.update(stage='failed',error=traceback.format_exc());save('STATUS.json',state);save('EXIT.json',dict(acquisition_complete=False,error=state['error']));raise
