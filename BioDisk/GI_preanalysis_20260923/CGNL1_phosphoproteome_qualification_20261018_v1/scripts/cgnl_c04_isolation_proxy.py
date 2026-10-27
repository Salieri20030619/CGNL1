from pathlib import Path
import json,xml.etree.ElementTree as E,struct,zlib,base64,datetime
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_search_pilot/run_v2_all_spectra';ns={'m':'http://psi.hupo.org/ms/mzml'}
def save(f,z):f.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
links=json.loads((W/'CGNL1_SCAN_LINKS.json').read_text());wanted={x['HCD']['scan'] for x in links}|{x['HCD']['refs'][0].rsplit('scan=',1)[1] for x in links}
save(W/'ISOLATION_PROXY_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),rule='UseactualmzMLisolationwindow andlinkedprecedingMS1;sumintensitynearselectedmz +/-10ppm / totalcentroidintensityinwindow',limitations='Selectedisotope proxy only,notoriginalSEQUESTisolation-specificity implementation orSPSinterferencecorrection',no_biological_comparison=True))
data={}
for _,x in E.iterparse(O/'C04_raw_pilot/phospho_fraction_1.mzML',events=['end']):
 if x.tag!='{'+ns['m']+'}spectrum':continue
 scan=x.attrib['id'].rsplit('scan=',1)[1]
 if scan in wanted:
  d={'windows':[{v.attrib['accession']:v.attrib['value'] for v in el.findall('m:cvParam',ns)} for el in x.findall('.//m:isolationWindow',ns)]}
  for b in x.findall('.//m:binaryDataArray',ns):
   acc={v.attrib['accession'] for v in b.findall('m:cvParam',ns)}
   if not acc&{'MS:1000514','MS:1000515'}:continue
   raw=base64.b64decode(b.findtext('m:binary','',ns));raw=zlib.decompress(raw) if 'MS:1000574' in acc else raw;fmt='d' if 'MS:1000523' in acc else 'f';size=8 if fmt=='d' else 4;d['mz' if 'MS:1000514' in acc else 'intensity']=list(struct.unpack('<'+fmt*(len(raw)//size),raw))
  data[scan]=d
 x.clear()
out=[]
for x in links:
 h=x['HCD'];m=data[h['scan']]['windows'][0];parent=h['refs'][0].rsplit('scan=',1)[1];c=float(m['MS:1000827']);lo=float(m['MS:1000828']);hi=float(m['MS:1000829']);d=data[parent];peaks=[[mz,it] for mz,it in zip(d['mz'],d['intensity']) if c-lo<=mz<=c+hi];total=sum(it for mz,it in peaks);signal=sum(it for mz,it in peaks if abs(mz-h['mz'])<=h['mz']*10e-6)
 out.append(dict(HCD_scan=h['scan'],MS1_scan=parent,isolation_center=c,lower_offset=lo,upper_offset=hi,selected_mz=h['mz'],window_centroids=peaks,total_intensity=total,selected_mz_intensity=signal,selected_peak_fraction=signal/total if total else None))
save(W/'ISOLATION_PROXY_RESULTS.json',out);print(json.dumps(out))
