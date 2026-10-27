from pathlib import Path
import json,sys,hashlib,datetime,re,subprocess,os,fcntl
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_totalprotein_pilot_20261026';P=O/'C04_search_pilot'
name=sys.argv[1];assert name in ['protein_set1_fraction_1.raw','protein_set2_fraction_1.raw']
d=W/Path(name).stem;out=d/'search_CID';out.mkdir(exist_ok=True);lock=(out/'search.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
assert not (out/'RUN.json').exists(), 'Existing run requires audit; do not restart'
accept=json.loads((W/(name+'.ACCEPTANCE.json')).read_text());assert accept['acquisition_pass'] and accept['MS2_activation_counts']=={'collision-induced dissociation':accept['counts']['2']}
mz=d/(Path(name).stem+'.mzML');assert sha(mz)==accept['mzML_sha256']
db=P/'mouse_plus_contaminants.fasta';assert sha(db)=='2c649a31e1afccced1a1ad0ce0a8e13676ceac44e918e3a540df9be720ddbcb7'
template=P/'run_v3_CID/comet.params'
frozen=json.loads((O/'C04_full_fraction_qualification/JOINT_IDENTIFICATION_SPEC_BEFORE_RESULTS.json').read_text())
assert sha(template)==frozen['CID_parameter_template_sha256']
text=template.read_text();params,n=re.subn(r'^variable_mod02\s*=.*$','variable_mod02 = 0.0 X 0 3 -1 0 0 0.0',text,flags=re.M);assert n==1
params,n=re.subn(r'^print_ascorepro_score\s*=.*$','print_ascorepro_score = 0',params,flags=re.M);assert n==1
active={line.split('=',1)[0].strip():line.split('=',1)[1].split('#',1)[0].strip() for line in params.splitlines() if '=' in line and not line.lstrip().startswith('#')}
assert active['activation_method']=='CID' and active['num_threads']=='4'
assert float(active['variable_mod02'].split()[0])==0
assert all(abs(float(v.split()[0])-79.966331)>1e-6 for k,v in active.items() if re.fullmatch(r'variable_mod\d+',k))
(out/'comet.params').write_text(params,encoding='utf-8')
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),file=name,mzML_sha256=accept['mzML_sha256'],database_sha256=sha(db),engine_sha256=sha(P/'comet.linux.exe'),template_sha256=sha(template),parameters_sha256=sha(out/'comet.params'),changes=['variable phosphorylation disabled for total proteome','phosphosite localization scoring disabled'],fixed_modifications='TMT Nterm/K and carbamidomethyl C unchanged',variable_modification='M oxidation only',mode='CID MS2',threads=4,scope='single fixed fraction pilot; not whole protein cohort',planned_PSM_QA='top1 PSM per scan/charge; Evalue ties;(D+1)/T reverse cumulative minimum; pilot diagnostic only',protein_FDR_not_established=True,CGNL1_absence_not_inferred_from_pilot=True,disease_tests=0)
save('SPEC_BEFORE_SEARCH.json',spec)
mem={line.split(':')[0]:line.split(':')[1].strip() for line in Path('/proc/meminfo').read_text().splitlines()};assert int(mem['MemAvailable'].split()[0])>16*1024**2
link=out/'input.mzML';assert not link.exists();link.symlink_to(mz)
cmd=[str(P/'comet.linux.exe'),'-P'+str(out/'comet.params'),'-N'+str(out/'psms'),str(link)]
with (out/'search.log').open('w',encoding='utf-8') as log:
 child=subprocess.Popen(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT);save('RUN.json',dict(pid=child.pid,parent_pid=os.getpid(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
files=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in out.glob('psms*') if p.is_file()]
save('SEARCH_EXIT.json',dict(exit=code,files=files,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pilot_QA_pending=True,whole_project_complete=False));assert code==0 and (out/'psms.txt').exists() and (out/'psms.pep.xml').exists();print(json.dumps(dict(file=name,exit=code,outputs=len(files),pilot_QA_pending=True)))
