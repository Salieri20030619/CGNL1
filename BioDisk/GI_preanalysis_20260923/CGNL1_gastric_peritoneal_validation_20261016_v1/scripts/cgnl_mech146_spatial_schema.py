from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,hashlib,zipfile,gzip,h5py,pandas as pd,numpy as np
from scipy import sparse,io
from threadpoolctl import threadpool_limits
from PIL import Image
O=Path('C:/codex_test/output/CGNL1_spatial_public_20261016');p=O/'Human_spatial.zip';r=next(x for x in json.loads((O/'DOWNLOAD_RECEIPTS_140.json').read_text(encoding='utf-8')) if x['name']==p.name);assert r['success'] and p.stat().st_size==r['size'] and hashlib.md5(p.read_bytes()).hexdigest()==r['computed_md5'];E=O/'extracted_spatial';E.mkdir(exist_ok=True)
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 for m in z.infolist():
  target=(E/m.filename).resolve();assert target.is_relative_to(E.resolve())
  if m.is_dir():target.mkdir(parents=True,exist_ok=True)
  else:
   target.parent.mkdir(parents=True,exist_ok=True);b=z.read(m)
   if target.exists():assert target.read_bytes()==b
   else:target.write_bytes(b)
D=E/'data';F=D/'filtered_feature_bc_matrix';out={'archive_md5':r['computed_md5'],'duplicate_files':{}}
for name in ['features.tsv','barcodes.tsv','matrix.mtx']:
 a=(F/name).read_bytes();b=gzip.decompress((F/(name+'.gz')).read_bytes());out['duplicate_files'][name]={'plain_bytes':len(a),'gzip_decoded_bytes':len(b),'byte_identical':a==b,'plain_sha256':hashlib.sha256(a).hexdigest(),'gzip_decoded_sha256':hashlib.sha256(b).hexdigest()}
with h5py.File(D/'filtered_feature_bc_matrix.h5','r') as h:
 g=h['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]);bar=np.array([x.decode('utf-8') for x in g['barcodes'][:]]);ids=np.array([x.decode('utf-8') for x in g['features/id'][:]]);genes=np.array([x.decode('utf-8') for x in g['features/name'][:]]);out['h5']={'shape':list(map(int,X.shape)),'dtype':str(X.dtype),'nonzero':int(X.nnz),'total_UMI':float(X.sum()),'unique_barcodes':len(set(bar))==len(bar),'features':list(g['features'].keys()),'target_rows':np.where(genes=='CGNL1')[0].tolist(),'target_ids':ids[genes=='CGNL1'].tolist()}
assert len(bar)==len(set(bar));assert X.data.min()>=0 and np.equal(X.data,np.floor(X.data)).all();assert sum(genes=='CGNL1')==1 and ids[genes=='CGNL1'][0]=='ENSG00000128849'
positions=pd.read_csv(D/'spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']);assert positions.barcode.is_unique;assert set(bar)<=set(positions.barcode);pos=positions.set_index('barcode').loc[bar];sf=json.loads((D/'spatial/scalefactors_json.json').read_text());out['coordinates']={'rows':len(positions),'matrix_spots':len(bar),'in_tissue_in_matrix':int(pos.in_tissue.sum()),'all_matrix_spots_have_coordinates':True,'scalefactors':sf,'images':{}}
for n,scale in [('tissue_hires_image.png','tissue_hires_scalef'),('tissue_lowres_image.png','tissue_lowres_scalef')]:
 with Image.open(D/'spatial'/n) as image:
  width,height=image.size;x=pos.pxl_col_in_fullres.to_numpy()*sf[scale];y=pos.pxl_row_in_fullres.to_numpy()*sf[scale];within=(x>=0)&(x<width)&(y>=0)&(y<height);assert within.all();out['coordinates']['images'][n]={'size':[width,height],'all_spot_centers_inside_image':True}
with threadpool_limits(limits=2):M=io.mmread(F/'matrix.mtx').tocsc()
plain_bar=pd.read_csv(F/'barcodes.tsv',header=None)[0].to_numpy();ft=pd.read_csv(F/'features.tsv',sep='\t',header=None);out['matrix_consistency']={'plain_shape':list(M.shape),'plain_barcodes_match_h5':np.array_equal(plain_bar,bar),'plain_ids_match_h5':np.array_equal(ft[0].to_numpy(),ids)}
if M.shape==X.shape and out['matrix_consistency']['plain_barcodes_match_h5'] and out['matrix_consistency']['plain_ids_match_h5']:out['matrix_consistency']['plain_counts_equal_h5']=(M!=X).nnz==0
(O/'SPATIAL_SCHEMA_146.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
