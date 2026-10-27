from pathlib import Path
import sys,json,re,hashlib,subprocess,datetime,os
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');F=O/'C04_full_fraction_qualification';P=O/'C04_search_pilot'
def read(f):return json.loads(f.read_text())
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(f,z):f.write_text(json.dumps(z,indent=2)+'\n',encoding='utf-8');assert read(f)==z
name=sys.argv[1];spec=read(F/'JOINT_IDENTIFICATION_SPEC_BEFORE_RESULTS.json');assert name in spec['files'];assert re.fullmatch(r'phospho_(fraction_\d+|pY)\.raw',name)
assert sha(P/'mouse_plus_contaminants.fasta')==spec['database_sha256']
if name=='phospho_fraction_1.raw':
 for sub in ['run_v2_all_spectra','run_v3_CID']:assert read(P/sub/'SEARCH_EXIT.json')['exit']==0
 save(F/'REUSED_FIRST_FRACTION_SEARCH.json',dict(name=name,HCD=str(P/'run_v2_all_spectra'),CID=str(P/'run_v3_CID'),reused=True,cohort_q_not_computed=True));print('Firstfractionsearch reused');raise SystemExit
stem=Path(name).stem;C=F/'converted'/stem;accept=read(C/'ACCEPTANCE.json');assert accept['acquisition_pass'];mz=C/(stem+'.mzML');assert sha(mz)==accept['mzML_sha256']
records=[]
for mode,template,required in [('HCD',P/'run_v2_all_spectra/comet.params','beam-type collision-induced dissociation'),('CID',P/'run_v3_CID/comet.params','collision-induced dissociation')]:
 expected=spec[mode+'_parameter_template_sha256'];assert sha(template)==expected
 w=F/'searches'/stem/mode;w.mkdir(parents=True,exist_ok=True)
 if not accept['MS2_activation_counts'].get(required,0):
  z=dict(name=name,mode=mode,status='no_corresponding_acquisition_mode',not_gene_absence=True);save(w/'NOT_APPLICABLE.json',z);records.append(z);continue
 if (w/'SEARCH_EXIT.json').exists() and read(w/'SEARCH_EXIT.json')['exit']==0:
  assert (w/'psms.txt').exists() and (w/'psms.pep.xml').exists();records.append(dict(mode=mode,reused_completed=True));continue
 if (w/'RUN.json').exists():
  old=read(w/'RUN.json');proc=Path('/proc')/str(old['pid'])
  if proc.exists() and 'comet.linux.exe' in (proc/'cmdline').read_bytes().decode().replace('\x00',' '):raise RuntimeError('Existing search process live;adopt it,do not duplicate')
  raise RuntimeError('Prior incomplete search requires explicit diagnosis/versioned recovery')
 (w/'comet.params').write_text(template.read_text(),encoding='utf-8');link=w/'input.mzML'
 if not link.exists():link.symlink_to(mz)
 cmd=[str(P/'comet.linux.exe'),'-P'+str(w/'comet.params'),'-N'+str(w/'psms'),str(link)]
 save(w/'INPUT_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),name=name,mode=mode,params_sha256=expected,mzML_sha256=accept['mzML_sha256'],database_sha256=spec['database_sha256'],family=('pY_' if name=='phospho_pY.raw' else 'ST_')+mode))
 available_kb=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))
 if available_kb<16*1024*1024:raise RuntimeError('Insufficient available memory for bounded search; raw/conversion preserved')
 with (w/'search.log').open('w',encoding='utf-8') as log:
  p=subprocess.Popen(cmd,cwd=w,stdout=log,stderr=subprocess.STDOUT);save(w/'RUN.json',dict(pid=p.pid,parent_pid=os.getpid(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=p.wait()
 z=dict(name=name,mode=mode,exit=code,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in w.glob('psms*') if f.is_file()]);save(w/'SEARCH_EXIT.json',z);assert code==0;records.append(z)
save(F/'searches'/stem/'FILE_SEARCH_COMPLETE.json',dict(name=name,records=records,cohort_q_not_computed=True,no_biological_inference=True));print(json.dumps(dict(name=name,file_search_complete=True)))
