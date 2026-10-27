from pathlib import Path
import json, hashlib, datetime, xml.etree.ElementTree as E, base64, zlib, struct

P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1')
F=P/'C04_full_fraction_qualification'; A=F/'joint_identification'; W=F/'joint_spectral_qualification_20261026'
W.mkdir(exist_ok=True)
ns={'m':'http://psi.hupo.org/ms/mzml'}; tag='{'+ns['m']+'}spectrum'
targets=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145]
def save(name,x):
 (W/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for b in iter(lambda:stream.read(8388608),b''):h.update(b)
 return h.hexdigest()
annotations=json.loads((A/'CGNL1_SITE_ANNOTATIONS.json').read_text())['candidates']
passed=[x for x in annotations if x['PSM_dual_q_pass']]
assert len(passed)==6
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),annotation_sha256=sha(A/'CGNL1_SITE_ANNOTATIONS.json'),scope='all six joint dual-q passing CGNL1 PSMs',pair_rule='same nonempty MS1 reference, charge, m/z within 10 ppm and RT within 2 seconds; unique opposite activation spectrum',MS3_rule='explicit CID precursor spectrumRef; no adjacency inference',reporters=targets,reporter_tolerance_ppm=20,imputation=False,disease_tests=0,replication_of_existing_pilot_rules=True)
assert not (W/'SPEC_BEFORE_EXTRACTION.json').exists(), 'Existing run requires audit'
save('SPEC_BEFORE_EXTRACTION.json',spec)
results=[]; manifests=[]
for name in sorted({x['file'] for x in passed}):
 stem=Path(name).stem
 path=P/'C04_raw_pilot/phospho_fraction_1.mzML' if name=='phospho_fraction_1.raw' else F/'converted'/stem/(stem+'.mzML')
 receipt=P/'C04_raw_pilot/ACQUISITION_ACCEPTANCE.json' if name=='phospho_fraction_1.raw' else F/'converted'/stem/'ACCEPTANCE.json'
 digest=sha(path); expected=json.loads(receipt.read_text()).get('mzML_sha256')
 if name=='phospho_fraction_1.raw':expected='adab96e2d5db2be34eedc1dc6529c1799a65a5b72314693a7b65b44c1b609449'
 assert digest==expected
 manifests.append(dict(file=name,mzML_path=str(path),sha256=digest))
 meta=[]
 for _,x in E.iterparse(path,events=['end']):
  if x.tag!=tag:continue
  level=int(next(v.attrib['value'] for v in x.findall('m:cvParam',ns) if v.attrib['accession']=='MS:1000511'))
  if level in [2,3]:
   cv=x.findall('.//m:cvParam',ns);rt=next(v for v in cv if v.attrib['accession']=='MS:1000016')
   ion=x.find('m:precursorList/m:precursor/m:selectedIonList/m:selectedIon',ns)
   vals={} if ion is None else {v.attrib['accession']:v.attrib.get('value') for v in ion.findall('m:cvParam',ns)}
   modes=sorted(set(v.attrib['name'] for v in x.findall('.//m:activation/m:cvParam',ns) if 'dissociation' in v.attrib.get('name','')))
   meta.append(dict(id=x.attrib['id'],scan=int(x.attrib['id'].rsplit('scan=',1)[1]),level=level,rt=float(rt.attrib['value'])*(60 if rt.attrib.get('unitAccession')=='UO:0000031' else 1),mz=float(vals.get('MS:1000744','0')),charge=vals.get('MS:1000041'),refs=sorted({v.attrib['spectrumRef'] for v in x.findall('m:precursorList/m:precursor',ns) if 'spectrumRef' in v.attrib}),mode='HCD' if modes==['beam-type collision-induced dissociation'] else 'CID' if modes==['collision-induced dissociation'] else str(modes)))
  x.clear()
 byscan={x['scan']:x for x in meta};groups={}
 for seed in [x for x in passed if x['file']==name]:
  s=byscan[seed['scan']];assert s['mode']==seed['mode'] and int(s['charge'])==seed['charge']
  opposite='CID' if s['mode']=='HCD' else 'HCD'
  mates=[x for x in meta if x['level']==2 and x['mode']==opposite and s['refs'] and x['refs']==s['refs'] and x['charge']==s['charge'] and abs(x['mz']-s['mz'])<=s['mz']*10e-6 and abs(x['rt']-s['rt'])<=2]
  h=s if s['mode']=='HCD' else mates[0] if len(mates)==1 else None
  c=s if s['mode']=='CID' else mates[0] if len(mates)==1 else None
  key=(h['scan'] if h else None,c['scan'] if c else None)
  if key in groups:continue
  linked=[x for x in meta if x['level']==3 and c and c['id'] in x['refs']]
  anns=[x for x in annotations if x['file']==name and x['scan'] in key]
  groups[key]=dict(file=name,HCD=h,CID=c,opposite_match_count=len(mates),MS3=linked,PSM_annotations=anns,paired_peptide_agreement=len(anns)==2 and len({x['peptide'] for x in anns})==1,paired_proposed_site_agreement=len(anns)==2 and anns[0]['engine_proposed_canonical_sites']==anns[1]['engine_proposed_canonical_sites'])
 selected={x['scan'] for g in groups.values() for x in g['MS3']};spectra={}
 for _,node in E.iterparse(path,events=['end']):
  if node.tag!=tag:continue
  scan=int(node.attrib['id'].rsplit('scan=',1)[1])
  if scan in selected:
   arrays={}
   for b in node.findall('.//m:binaryDataArray',ns):
    cvs={x.attrib['accession'] for x in b.findall('m:cvParam',ns)}
    if not cvs&{'MS:1000514','MS:1000515'}:continue
    assert ('MS:1000523' in cvs) != ('MS:1000521' in cvs)
    raw=base64.b64decode(b.findtext('m:binary','',ns),validate=True);raw=zlib.decompress(raw) if 'MS:1000574' in cvs else raw
    fmt='d' if 'MS:1000523' in cvs else 'f';size=8 if fmt=='d' else 4;assert len(raw)%size==0
    arrays['mz' if 'MS:1000514' in cvs else 'intensity']=struct.unpack('<'+fmt*(len(raw)//size),raw)
   assert len(arrays['mz'])==len(arrays['intensity'])
   peaks=[]
   for target in targets:
    matches=sorted((abs(m-target),m,i) for m,i in zip(arrays['mz'],arrays['intensity']) if abs(m-target)<=target*20e-6)
    peaks.append(dict(target_mz=target,matches=len(matches),observed_mz=matches[0][1] if matches else None,intensity=matches[0][2] if matches else None))
   spectra[scan]=peaks
  node.clear()
 assert set(spectra)==selected
 for g in groups.values():
  for x in g['MS3']:x['reporter_peaks']=spectra[x['scan']]
 results.extend(groups.values())
save('RESULTS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputs=manifests,groups=results,biological_tests=0,original_SNR_and_PIF_equivalence_not_established=True,total_protein_gate_pending=True,whole_project_complete=False))
print(json.dumps([dict(file=g['file'],HCD=g['HCD']['scan'] if g['HCD'] else None,CID=g['CID']['scan'] if g['CID'] else None,peptide_agreement=g['paired_peptide_agreement'],site_assignment_agreement=g['paired_proposed_site_agreement'],MS3=[dict(scan=x['scan'],positive_unique_reporters=sum(r['matches']==1 and r['intensity']>0 for r in x['reporter_peaks'])) for x in g['MS3']]) for g in results]))
