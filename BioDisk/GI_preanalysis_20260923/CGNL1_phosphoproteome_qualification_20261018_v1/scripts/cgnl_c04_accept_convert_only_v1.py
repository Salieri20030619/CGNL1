from pathlib import Path
import sys,re,json,hashlib,subprocess,datetime,xml.etree.ElementTree as E,collections
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_full_fraction_qualification';R=W/'raw';R.mkdir(exist_ok=True)
name,digest=sys.argv[1:];assert re.fullmatch(r'phospho_(fraction_\d+|pY)\.raw',name)
spec=json.loads((W/'ACQUISITION_SPEC_BEFORE_RESULTS.json').read_text());entry=next(x for x in spec['files'] if x['name']==name);f=R/name;upload=R/(name+'.uploading');out=W/'converted'/f.stem;out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,z):p.write_text(json.dumps(z,indent=2)+'\n');assert json.loads(p.read_text())==z
if (out/'ACCEPTANCE.json').exists():
 old=json.loads((out/'ACCEPTANCE.json').read_text());assert old['raw_sha256']==digest and sha(f)==digest;print(json.dumps(dict(reused=True,**old)));raise SystemExit
if upload.exists():
 assert upload.stat().st_size==entry['bytes'] and sha(upload)==digest
 assert not f.exists(),'Existingunacceptedfinalrequiresinspection'
 upload.rename(f)
assert f.stat().st_size==entry['bytes'] and sha(f)==digest
save(out/'INPUT_ACCEPTANCE.json',dict(name=name,bytes=f.stat().st_size,sha256=digest,utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
binary=O/'software/ThermoRawFileParser_2_dev/ThermoRawFileParser'
with (out/'conversion.log').open('w') as log:
 p=subprocess.Popen([str(binary),'-i='+str(f),'-o='+str(out),'-f=2','-m=0'],stdout=log,stderr=subprocess.STDOUT)
 save(out/'RUN.json',dict(pid=p.pid,utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=p.wait()
save(out/'CONVERSION_EXIT.json',dict(exit=code));assert code==0
mz=out/(f.stem+'.mzML');ns={'m':'http://psi.hupo.org/ms/mzml'};counts=collections.Counter();modes=collections.Counter();ids=set();badrefs=0;ms3refs=0
for _,x in E.iterparse(mz,events=['end']):
 if x.tag!='{'+ns['m']+'}spectrum':continue
 level=int(next(v.attrib['value'] for v in x.findall('m:cvParam',ns) if v.attrib['accession']=='MS:1000511'));counts[level]+=1
 if level==2:
  mode=';'.join(sorted(set(v.attrib['name'] for v in x.findall('.//m:activation/m:cvParam',ns) if 'dissociation' in v.attrib.get('name',''))));modes[mode]+=1
 if level==3:
  refs=[v.attrib.get('spectrumRef') for v in x.findall('.//m:precursor',ns)];ms3refs+=len(refs);badrefs+=sum(r not in ids for r in refs)
 ids.add(x.attrib['id']);x.clear()
log=(out/'conversion.log').read_text();warnings=[x for x in log.splitlines() if 'WARN' in x]
z=dict(name=name,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),raw_bytes=f.stat().st_size,raw_sha256=digest,mzML_bytes=mz.stat().st_size,mzML_sha256=sha(mz),spectrum_counts={str(k):v for k,v in counts.items()},MS2_activation_counts=dict(modes),MS3_reference_count=ms3refs,unresolved_MS3_references=badrefs,conversion_warnings=warnings,conversion_exit=code,acquisition_pass=counts[2]>0 and badrefs==0,MS3_quantification_present=counts[3]>0,CGNL1_measured=None,biological_analysis_complete=False)
save(out/'ACCEPTANCE.json',z);print(json.dumps(z))
