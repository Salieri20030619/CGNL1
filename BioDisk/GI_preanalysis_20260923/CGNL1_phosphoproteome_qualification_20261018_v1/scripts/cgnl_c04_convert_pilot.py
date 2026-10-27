from pathlib import Path
import json,hashlib,subprocess,datetime,xml.etree.ElementTree as E,collections,base64,zlib,struct
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');D=O/'C04_sources';W=O/'C04_raw_pilot';W.mkdir(exist_ok=True)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,x):f.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==x
raw=D/'phospho_fraction_1.raw';receipt=json.loads((D/'C04_FTPS_RAW_RECEIPT.json').read_text());assert raw.stat().st_size==221702430 and sha(raw)==receipt['sha256']
binary=O/'software/ThermoRawFileParser_2_dev/ThermoRawFileParser'
save(W/'CONVERSION_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),raw_sha256=sha(raw),binary_sha256=sha(binary),method='IndexedmzML includingallMSlevels,defaultvendorpeakpicking',checks=['RAWreadexit','MS1/2/3counts','MS3precursorreferences','MS3reportermz +/-20ppm at8assignedTMTchannels'],reporter_mz=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145],no_peptide_identification=True))
with (W/'conversion.log').open('w') as log:
 r=subprocess.run([str(binary),'-i='+str(raw),'-o='+str(W),'-f=2','-m=0'],stdout=log,stderr=subprocess.STDOUT)
save(W/'CONVERSION_EXIT.json',dict(exit=r.returncode,utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
assert r.returncode==0,(W/'conversion.log').read_text()[-1000:]
f=W/'phospho_fraction_1.mzML';assert f.exists()
ns={'m':'http://psi.hupo.org/ms/mzml'};counts=collections.Counter();ids=set();ms3refs=0;missingrefs=[];withreporter=0;channelcounts=[0]*8;mz_targets=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145]
def array(node):
 cvs={x.attrib['accession'] for x in node.findall('m:cvParam',ns)};data=base64.b64decode(node.findtext('m:binary','',ns))
 if 'MS:1000574' in cvs:data=zlib.decompress(data)
 fmt='d' if 'MS:1000523' in cvs else 'f';size=8 if fmt=='d' else 4
 assert len(data)%size==0;return struct.unpack('<'+fmt*(len(data)//size),data)
for event,node in E.iterparse(f,events=['end']):
 if node.tag!='{'+ns['m']+'}spectrum':continue
 cv={x.attrib['accession']:x.attrib.get('value') for x in node.findall('m:cvParam',ns)};level=int(cv['MS:1000511']);counts[level]+=1
 if level==3:
  refs=[x.attrib.get('spectrumRef') for x in node.findall('.//m:precursor',ns)];ms3refs+=len(refs)
  missingrefs.extend(x for x in refs if x not in ids)
  arrays={}
  for b in node.findall('.//m:binaryDataArray',ns):
   acc={x.attrib['accession'] for x in b.findall('m:cvParam',ns)}
   if 'MS:1000514' in acc:arrays['mz']=array(b)
   if 'MS:1000515' in acc:arrays['intensity']=array(b)
  if 'mz' in arrays and 'intensity' in arrays:
   hit=False
   for j,target in enumerate(mz_targets):
    yes=any(abs(m-target)<=target*20e-6 and intensity>0 for m,intensity in zip(arrays['mz'],arrays['intensity']))
    channelcounts[j]+=int(yes);hit|=yes
   withreporter+=int(hit)
 ids.add(node.attrib['id']);node.clear()
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),raw_sha256=sha(raw),mzML_sha256=sha(f),mzML_bytes=f.stat().st_size,conversion_exit=0,spectrum_counts={str(k):v for k,v in counts.items()},ms3_precursor_references=ms3refs,unresolved_precursor_refs=len(missingrefs),ms3_with_any_assigned_reporter_peak=withreporter,reporter_peak_scan_counts=channelcounts,reporter_mz=mz_targets,CGNL1_identified=None,biological_analysis_complete=False,search_localization_not_run=True)
save(W/'ACQUISITION_QA.json',z);print(json.dumps(z),flush=True)
