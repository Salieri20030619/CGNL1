import pathlib,sys,tarfile,json,numpy as np,pandas as pd
from scipy.spatial import cKDTree
O=pathlib.Path(sys.argv[1]);out=[]
for name in ['IX_56689','M-ST-11']:
 with tarfile.open(O/'raw/liver'/(name+'.tar.gz')) as t:
  names=t.getnames();sf=json.load(t.extractfile([x for x in names if x.endswith('scalefactors_json.json')][0]));p=[x for x in names if x.endswith('tissue_positions.csv')];p=p or [x for x in names if x.endswith('tissue_positions_list.csv')];a=pd.read_csv(t.extractfile(p[0]),header=0 if p[0].endswith('tissue_positions.csv') else None)
  if 'barcode' not in a.columns:a.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
  xy=a[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);out.append({'sample':name,'scalefactors':sf,'pixel_pitch':pitch,'pitch_um_from55':pitch*55/sf['spot_diameter_fullres'],'implied_diameter_using100pitch':sf['spot_diameter_fullres']*100/pitch,'rows':a.array_row.max(),'cols':a.array_col.max(),'spacing_from_grid':float(np.median(cKDTree(np.c_[a.array_col*50,a.array_row*50*np.sqrt(3)]).query(np.c_[a.array_col*50,a.array_row*50*np.sqrt(3)],k=2)[0][:,1]))})
(O/'results/coordinate_calibration_audit.json').write_text(json.dumps(out,indent=2,default=lambda x: x.item()));print(json.dumps(out,indent=2,default=lambda x: x.item()))
