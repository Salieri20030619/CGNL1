from pathlib import Path
import json,gzip,hashlib,zipfile,io,re,datetime
import numpy as np,pandas as pd,h5py
from scipy.spatial import cKDTree
from scipy import sparse,io as sio
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parents[1];B=R.parent;A=B/'expanded_screen_20260925';D=R/'03_spatial_geometry'
def canon(X,genes,bc):
 assert len(set(bc))==len(bc)
 unique,inverse=np.unique(genes,return_inverse=True)
 collapse=sparse.csr_matrix((np.ones(len(genes)),(inverse,np.arange(len(genes)))),shape=(len(unique),len(genes)))
 X=(collapse@X).tocsr();genes=unique
 gi=np.argsort(genes);ci=np.argsort(bc);X=X[gi][:,ci].tocsr();X.eliminate_zeros();X.sort_indices();assert np.allclose(X.data,np.round(X.data))
 h=hashlib.sha256()
 for arr in [np.array(X.shape,dtype='<i8'),X.indptr.astype('<i8'),X.indices.astype('<i8'),X.data.astype('<i8')]:h.update(arr.tobytes())
 h.update('\n'.join(genes[gi]).encode());h.update('\n'.join(bc[ci]).encode());return h.hexdigest()
def h5read(f):
 with h5py.File(f,'r') as h:
  m=h['matrix'];X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));genes=m['features/name'].asstr()[:];bc=m['barcodes'].asstr()[:]
 return X,genes,bc
rows=[];fingerprints=[];patient=[]
for ds in ['GSE206552','GSE254364','GSE217414']:
 p=A/'spatial_data'/ds
 for s in json.loads((A/'metadata'/f'{ds}_samples.json').read_text()):
  gsm=s['gsm'];candidates=[x for x in p.glob(gsm+'*') if x.is_file()]+[x for x in (p/(gsm+'_unpacked')).rglob('*') if x.is_file()]
  def find(term):
   z=[x for x in candidates if term in x.name];assert len(z)==1,(gsm,term,z);return z[0]
  sfp=find('scalefactors');sf=json.loads(gzip.decompress(sfp.read_bytes()) if sfp.suffix=='.gz' else sfp.read_bytes());diam=float(sf['spot_diameter_fullres']);assert diam>0
  pos=pd.read_csv(find('tissue_positions'),header=None)
  if pos.iloc[0,0]=='barcode':pos=pos.iloc[1:]
  pos.columns=['barcode','in_tissue','array_row','array_col','y','x'];pos=pos.set_index('barcode').apply(pd.to_numeric);assert pos.index.is_unique
  xy=pos[['x','y']].to_numpy();ratio=55/diam;nn=cKDTree(xy).query(xy,k=2)[0][:,1]*ratio;med=float(np.median(nn));pitch_ok=80<=med<=120
  pos['x_um']=pos.x*ratio;pos['y_um']=pos.y*ratio;pos.to_csv(D/f'{gsm}_physical_coordinates.tsv.gz',sep='\t')
  h5=[x for x in candidates if x.name.endswith('_filtered_feature_bc_matrix.h5')]
  if h5:X,genes,bc=h5read(h5[0])
  else:
   genes=pd.read_csv(find('features.tsv'),sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(find('barcodes.tsv'),sep='\t',header=None)[0].to_numpy(dtype=str)
   with gzip.open(find('matrix.mtx'),'rb') as f,threadpool_limits(limits=2):X=sio.mmread(f).tocsc()
  assert set(bc)<=set(pos.index);tissue=pos.loc[bc].in_tissue.eq(1)
  fingerprint=canon(X,genes,bc);fingerprints.append(dict(dataset=ds,sample=gsm,sha256_canonical_matrix=fingerprint,source=str(h5[0]) if h5 else str(find('matrix.mtx'))))
  for g in ['CGNL1','CLDN5','ESAM','TJP1','OCLN']:assert g in genes,(gsm,g)
  rows.append(dict(dataset=ds,sample=gsm,spots=len(bc),tissue_spots=int(tissue.sum()),um_per_fullres_pixel=ratio,nearest_pitch_um=med,pitch_check=pitch_ok,extent_x_um=float(np.ptp(pos.loc[bc].x_um)),extent_y_um=float(np.ptp(pos.loc[bc].y_um)),panel_coverage='5/5 including target',pathology_region_status='pending_target_blind_region_review',fixed_distance_endpoint_status='not_yet_qualified'))
  match=re.search(r'Patient\s+(\d+)',s['title']);patient.append(dict(dataset=ds,sample=gsm,patient_id=ds+':'+match.group(1) if match else '',patient_evidence='author_GEO_title' if match else 'unresolved',title=s['title'],characteristics=json.dumps(s.get('characteristics',[])),cross_study_independence='pending_original_source_review'))
  print('GEOMETRY',gsm,round(med,3),flush=True)
p=B/'spatial_v3/background_validation/run_20260924T052010Z/inputs/NicheDE'
for f in sorted(p.glob('*.zip')):
 with zipfile.ZipFile(f) as z:
  names=[n for n in z.namelist() if not n.startswith('__MACOSX/') and n.endswith('/filtered_feature_bc_matrix.h5')]
  assert len(names)==1
  X,genes,bc=h5read(io.BytesIO(z.read(names[0])));fingerprints.append(dict(dataset='NicheDE_existing',sample=f.stem,sha256_canonical_matrix=canon(X,genes,bc),source=str(f)+'::'+names[0]))
fp=pd.DataFrame(fingerprints);fp.to_csv(D/'matrix_fingerprints.tsv',sep='\t',index=False)
duplicates=fp[fp.duplicated('sha256_canonical_matrix',keep=False)];duplicates.to_csv(D/'exact_duplicate_matrices.tsv',sep='\t',index=False)
pd.DataFrame(rows).to_csv(D/'geometry_qualification.tsv',sep='\t',index=False);pd.DataFrame(patient).to_csv(D/'patient_provenance.tsv',sep='\t',index=False)
summary={'exit_code':0,'status':'measurement_geometry_complete_region_eligibility_pending','sections':len(rows),'pitch_checks_passed':sum(x['pitch_check'] for x in rows),'matrices_compared':len(fp),'exact_duplicate_matrix_records':len(duplicates),'GEO_explicit_patient_labels':sum(bool(x['patient_id']) for x in patient),'patient_independence':'not established by matrix fingerprint; original-source review pending','target_tests_run':False,'finished':datetime.datetime.now().astimezone().isoformat()}
(D/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
