from pathlib import Path
import json,xml.etree.ElementTree as E,base64,zlib,struct,bisect,datetime,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=P/'C04_full_fraction_qualification/joint_spectral_qualification_20261026'
ns={'m':'http://psi.hupo.org/ms/mzml'};tag='{'+ns['m']+'}spectrum'
cvmap={'MS:1000514':'mz','MS:1000515':'intensity','MS:1002743':'noise_mz','MS:1002744':'noise','MS:1002745':'baseline'}
targets=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145]
def save(n,x):(W/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def interp(xs,ys,x):
 i=bisect.bisect_left(xs,x)
 if i<len(xs) and xs[i]==x:return ys[i]
 if i==0 or i==len(xs):return None
 return ys[i-1]+(ys[i]-ys[i-1])*(x-xs[i-1])/(xs[i]-xs[i-1])
assert not (W/'NOISE_SPEC.json').exists()
save('NOISE_SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope={'fraction1':[26697,26749],'fraction11':[19264,19300]},reporter_ppm=20,noise_mapping='exact sampled noise m/z or linear interpolation between bracketing noise samples; no extrapolation',metrics=['I/noise','(I-baseline)/noise'],original_SNR_equivalence=False,disease_tests=0,missing_reporters_not_imputed=True))
outputs=[];raw=[]
for name,path,scans in [('fraction1',P/'C04_raw_pilot/noise_MS3/phospho_fraction_1.mzML',{26697,26749}),('fraction11',W/'noise_fraction11/phospho_fraction_11.mzML',{19264,19300})]:
 seen=set()
 for _,x in E.iterparse(path,events=['end']):
  if x.tag!=tag:continue
  scan=int(x.attrib['id'].rsplit('scan=',1)[1])
  if scan in scans:
   arrays={}
   for b in x.findall('.//m:binaryDataArray',ns):
    cvs={v.attrib['accession'] for v in b.findall('m:cvParam',ns)};keys=cvs&set(cvmap)
    if not keys:continue
    assert len(keys)==1
    assert ('MS:1000523' in cvs)!=('MS:1000521' in cvs)
    data=base64.b64decode(b.findtext('m:binary','',ns),validate=True);data=zlib.decompress(data) if 'MS:1000574' in cvs else data;fmt='d' if 'MS:1000523' in cvs else 'f';size=8 if fmt=='d' else 4;assert len(data)%size==0
    arrays[cvmap[next(iter(keys))]]=list(struct.unpack('<'+fmt*(len(data)//size),data))
   assert set(arrays)==set(cvmap.values())
   assert len(arrays['mz'])==len(arrays['intensity']) and len(arrays['noise_mz'])==len(arrays['noise'])==len(arrays['baseline'])
   assert all(a<b for a,b in zip(arrays['noise_mz'],arrays['noise_mz'][1:]))
   reporters=[]
   for target in targets:
    matches=[i for i,m in enumerate(arrays['mz']) if abs(m-target)<=target*20e-6]
    row=dict(target_mz=target,matches=len(matches),I_over_noise=None,baseline_adjusted_proxy=None)
    if len(matches)==1:
     i=matches[0];m=arrays['mz'][i];v=arrays['intensity'][i];noise=interp(arrays['noise_mz'],arrays['noise'],m);baseline=interp(arrays['noise_mz'],arrays['baseline'],m)
     row.update(observed_mz=m,intensity=v,noise=noise,baseline=baseline)
     if noise is not None and noise>0:row.update(I_over_noise=v/noise,baseline_adjusted_proxy=(v-baseline)/noise)
    reporters.append(row)
   ratios=[r['I_over_noise'] for r in reporters if r['I_over_noise'] is not None]
   outputs.append(dict(file=name,scan=scan,reporters=reporters,available_reporter_ratios=len(ratios),sum_available_I_over_noise=sum(ratios),all_eight_evaluable=len(ratios)==8))
   raw.append(dict(file=name,scan=scan,arrays=arrays));seen.add(scan)
  x.clear()
 assert seen==scans
save('TARGET_NOISE_ARRAYS.json',raw);save('NOISE_RESULTS.json',dict(records=outputs,original_SNR_equivalence=False,new_biological_tests=0,whole_project_complete=False))
print(json.dumps([{k:r[k] for k in ['file','scan','available_reporter_ratios','sum_available_I_over_noise','all_eight_evaluable']} for r in outputs]))
