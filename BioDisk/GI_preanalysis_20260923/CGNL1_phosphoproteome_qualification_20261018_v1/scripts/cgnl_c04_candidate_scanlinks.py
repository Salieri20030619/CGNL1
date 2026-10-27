from pathlib import Path
import json,xml.etree.ElementTree as E,datetime
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_search_pilot/run_v2_all_spectra';ns={'m':'http://psi.hupo.org/ms/mzml'}
def save(f,z):f.write_text(json.dumps(z,indent=2)+'\n');assert json.loads(f.read_text())==z
save(W/'SCAN_LINK_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),rule='HCDcandidate to uniqueCID withsameMS1precursorreference,samecharge,mzwithin10ppm,RTwithin2sec;thenMS3explicitprecursorrefs',no_adjacency_only_matching=True,no_disease_comparison=True))
candidates=json.loads((W/'CGNL1_CANDIDATE_PSMS.json').read_text());meta=[]
for _,x in E.iterparse(O/'C04_raw_pilot/phospho_fraction_1.mzML',events=['end']):
 if x.tag!='{'+ns['m']+'}spectrum':continue
 level=int(next(v.attrib['value'] for v in x.findall('m:cvParam',ns) if v.attrib['accession']=='MS:1000511'))
 if level in [2,3]:
  cvs=x.findall('.//m:cvParam',ns);rt=next(v for v in cvs if v.attrib['accession']=='MS:1000016');rtsec=float(rt.attrib['value'])*(60 if rt.attrib.get('unitAccession')=='UO:0000031' else 1)
  ion=x.find('m:precursorList/m:precursor/m:selectedIonList/m:selectedIon',ns);vals={} if ion is None else {v.attrib['accession']:v.attrib.get('value') for v in ion.findall('m:cvParam',ns)}
  modes=sorted(set(v.attrib['name'] for v in x.findall('.//m:activation/m:cvParam',ns) if 'dissociation' in v.attrib.get('name','')))
  meta.append(dict(id=x.attrib['id'],scan=x.attrib['id'].rsplit('scan=',1)[1],level=level,rtsec=rtsec,mz=float(vals.get('MS:1000744','0')),charge=vals.get('MS:1000041'),refs=sorted(set(v.attrib.get('spectrumRef') for v in x.findall('m:precursorList/m:precursor',ns))),modes=modes))
 x.clear()
byid={x['scan']:x for x in meta};out=[]
for c in candidates:
 h=byid[c['scan']];cid=[x for x in meta if x['level']==2 and x['modes']==['collision-induced dissociation'] and x['refs']==h['refs'] and x['charge']==h['charge'] and abs(x['mz']-h['mz'])<=h['mz']*10e-6 and abs(x['rtsec']-h['rtsec'])<=2]
 ms3=[x for x in meta if x['level']==3 and len(cid)==1 and cid[0]['id'] in x['refs']]
 out.append(dict(HCD=h,CID_candidates=cid,MS3_linked=ms3,global_q=c['q_global'],phospho_q=c.get('q_phospho'),peptide=c['plain_peptide'],localization_score=c['ascore_sitescores']))
save(W/'CGNL1_SCAN_LINKS.json',out);print(json.dumps(out))
