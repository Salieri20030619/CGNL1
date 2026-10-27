from pathlib import Path
import json,datetime,hashlib,xml.etree.ElementTree as E,base64,zlib,struct
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=P/'C04_full_fraction_qualification/joint_spectral_qualification_20261026'
ns={'m':'http://psi.hupo.org/ms/mzml'};tag='{'+ns['m']+'}spectrum'
def save(n,x):(W/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
r=json.loads((W/'RESULTS.json').read_text());assert len(r['groups'])==4
assert not (W/'ISOLATION_SPEC.json').exists()
save('ISOLATION_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),input_linkage_sha256=hashlib.sha256((W/'RESULTS.json').read_bytes()).hexdigest(),scope='all four precursor pairs, no outcome-based selection',window='actual HCD isolation center and asymmetric offsets; explicit linked MS1',tolerance_ppm=10,metrics=['selected_mz_peak_fraction','isolation_center_peak_fraction','charge_spaced_lattice_fraction'],lattice='selected m/z + k * 1.0033548378/charge, k=-2..6, only within actual window; not validated isotope identity',original_PIF_equivalence=False,threshold_based_exclusion=False,disease_tests=0))
out=[]
for item in r['inputs']:
 groups=[g for g in r['groups'] if g['file']==item['file']];wanted=set()
 for g in groups:
  assert len(g['HCD']['refs'])==1
  wanted.add(g['HCD']['scan']);wanted.add(int(g['HCD']['refs'][0].rsplit('scan=',1)[1]))
 data={}
 for _,node in E.iterparse(item['mzML_path'],events=['end']):
  if node.tag!=tag:continue
  scan=int(node.attrib['id'].rsplit('scan=',1)[1])
  if scan in wanted:
   cvs={v.attrib['accession']:v.attrib.get('value') for v in node.findall('m:cvParam',ns)}
   d=dict(ms_level=int(cvs['MS:1000511']),windows=[{v.attrib['accession']:v.attrib['value'] for v in el.findall('m:cvParam',ns)} for el in node.findall('.//m:isolationWindow',ns)])
   if d['ms_level']==1:
    for b in node.findall('.//m:binaryDataArray',ns):
     acc={v.attrib['accession'] for v in b.findall('m:cvParam',ns)}
     if not acc&{'MS:1000514','MS:1000515'}:continue
     raw=base64.b64decode(b.findtext('m:binary','',ns),validate=True);raw=zlib.decompress(raw) if 'MS:1000574' in acc else raw
     assert ('MS:1000523' in acc)!=('MS:1000521' in acc)
     fmt='d' if 'MS:1000523' in acc else 'f';size=8 if fmt=='d' else 4;assert len(raw)%size==0
     d['mz' if 'MS:1000514' in acc else 'intensity']=struct.unpack('<'+fmt*(len(raw)//size),raw)
   data[scan]=d
  node.clear()
 assert set(data)==wanted
 for g in groups:
  h=g['HCD'];parent=int(h['refs'][0].rsplit('scan=',1)[1]);d=data[parent];assert d['ms_level']==1 and data[h['scan']]['ms_level']==2
  windows=data[h['scan']]['windows'];assert len(windows)==1
  window=windows[0];center=float(window['MS:1000827']);lo=float(window['MS:1000828']);hi=float(window['MS:1000829']);z=int(h['charge']);assert z>0
  peaks=[dict(mz=m,intensity=i) for m,i in zip(d['mz'],d['intensity']) if center-lo<=m<=center+hi];assert all(x['intensity']>=0 for x in peaks)
  lattice=[h['mz']+k*1.0033548378/z for k in range(-2,7)]
  for x in peaks:
   x['selected_match']=abs(x['mz']-h['mz'])<=h['mz']*10e-6
   x['center_match']=abs(x['mz']-center)<=center*10e-6
   x['lattice_match']=any(abs(x['mz']-t)<=t*10e-6 for t in lattice)
  total=sum(x['intensity'] for x in peaks)
  fractions={name:sum(x['intensity'] for x in peaks if x[key])/total if total else None for name,key in [('selected_mz_peak_fraction','selected_match'),('isolation_center_peak_fraction','center_match'),('charge_spaced_lattice_fraction','lattice_match')]}
  out.append(dict(file=g['file'],HCD=h['scan'],CID=g['CID']['scan'],MS1=parent,selected_mz=h['mz'],charge=z,isolation_center=center,lower_offset=lo,upper_offset=hi,window_centroids=peaks,total_intensity=total,**fractions))
save('ISOLATION_RESULTS.json',dict(records=out,original_PIF_equivalence=False,SPS_interference_not_established=True,new_biological_tests=0,whole_project_complete=False))
print(json.dumps([{k:v for k,v in x.items() if k!='window_centroids'} for x in out]))
