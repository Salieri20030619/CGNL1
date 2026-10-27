from pathlib import Path
import json,xml.etree.ElementTree as E,base64,zlib,struct,datetime
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_search_pilot/run_v2_all_spectra';ns={'m':'http://psi.hupo.org/ms/mzml'}
def save(f,z):f.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
targets=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145]
save(W/'SPECTRAL_QUALITY_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scans=[26694,26695,26697,26745,26747,26749],reporters=targets,matching='nearestcentroidwithin20ppm,recordmultiplematchcountsandintensity;no imputation',quality_only=True,no_group_comparison=True,limitations='Centroidintensity isnot originalSNR;isolationpurity/SPSinterference notestablished'))
selected={str(x) for x in [26694,26695,26697,26745,26747,26749]};records=[]
for _,node in E.iterparse(O/'C04_raw_pilot/phospho_fraction_1.mzML',events=['end']):
 if node.tag!='{'+ns['m']+'}spectrum':continue
 scan=node.attrib['id'].rsplit('scan=',1)[1]
 if scan in selected:
  cv={x.attrib['accession']:x.attrib.get('value') for x in node.findall('m:cvParam',ns)};arrays={};names=[]
  for b in node.findall('.//m:binaryDataArray',ns):
   cvs={x.attrib['accession'] for x in b.findall('m:cvParam',ns)};names.extend(x.attrib['name'] for x in b.findall('m:cvParam',ns))
   if not cvs&{'MS:1000514','MS:1000515'}:continue
   raw=base64.b64decode(b.findtext('m:binary','',ns));raw=zlib.decompress(raw) if 'MS:1000574' in cvs else raw;fmt='d' if 'MS:1000523' in cvs else 'f';size=8 if fmt=='d' else 4
   vals=list(struct.unpack('<'+fmt*(len(raw)//size),raw));arrays['mz' if 'MS:1000514' in cvs else 'intensity']=vals
  record=dict(scan=int(scan),ms_level=int(cv['MS:1000511']),array_annotations=names,arrays=arrays)
  if record['ms_level']==3:
   peaks=[]
   for target in targets:
    found=[(abs(m-target),m,i) for m,i in zip(arrays['mz'],arrays['intensity']) if abs(m-target)<=target*20e-6]
    found.sort();peaks.append(dict(target_mz=target,matches=len(found),observed_mz=found[0][1] if found else None,intensity=found[0][2] if found else None,ppm=(found[0][1]-target)/target*1e6 if found else None))
   record['reporter_peaks']=peaks
  records.append(record)
 node.clear()
assert len(records)==6
save(W/'CANDIDATE_SPECTRA.json',records)
summary=[dict(scan=x['scan'],ms_level=x['ms_level'],peaks=len(x['arrays']['mz']),positive_reporters=sum((r['intensity'] or 0)>0 for r in x.get('reporter_peaks',[])),max_abs_reporter_ppm=max((abs(r['ppm']) for r in x.get('reporter_peaks',[]) if r['ppm'] is not None),default=None),multiple_reporter_matches=sum(r['matches']>1 for r in x.get('reporter_peaks',[]))) for x in records]
save(W/'SPECTRAL_QUALITY_RESULTS.json',dict(scans=summary,original_SNR_available=False,isolation_specificity_verified=False,biological_tests=0));print(json.dumps(summary))
