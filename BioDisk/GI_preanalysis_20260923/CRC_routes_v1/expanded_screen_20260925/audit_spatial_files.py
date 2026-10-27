from pathlib import Path
import gzip,json,tarfile,numpy as np,pandas as pd,h5py
from scipy import sparse,io
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
A=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/expanded_screen_20260925');D=A/'spatial_data';O=A/'results';O.mkdir(exist_ok=True);rows=[];fig,axes=plt.subplots(2,5,figsize=(20,8));k=0
for g in ['GSE206552','GSE254364','GSE217414']:
 samples=json.loads((A/'metadata'/(g+'_samples.json')).read_text());R=D/g
 for s in samples:
  gsm=s['gsm'];files=list(R.glob(gsm+'*'));T=R/(gsm+'_unpacked');T.mkdir(exist_ok=True)
  for p in files:
   if p.name.endswith('_spatial.tar.gz'):
    with tarfile.open(p) as f:
     for m in f.getmembers():assert (T/m.name).resolve().is_relative_to(T.resolve())
     f.extractall(T,filter='data')
  candidates=[p for p in files if p.is_file()]+list(T.rglob('*'));candidates=[p for p in candidates if p.is_file()]
  def find(pattern):
   z=[p for p in candidates if pattern in p.name];assert len(z)==1,(gsm,pattern,[str(x) for x in z]);return z[0]
  h=[p for p in files if p.name.endswith('_filtered_feature_bc_matrix.h5')]
  if h:
   with h5py.File(h[0]) as f:
    m=f['matrix'];X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));sy=m['features']['name'][:].astype(str);bc=m['barcodes'][:].astype(str)
  else:
   sy=pd.read_csv(find('features.tsv'),sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(find('barcodes.tsv'),sep='\t',header=None)[0].to_numpy(dtype=str)
   with gzip.open(find('matrix.mtx'),'rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsc()
  assert X.shape==(len(sy),len(bc)) and pd.Index(bc).is_unique;assert (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
  pp=find('tissue_positions');pos=pd.read_csv(pp,header=None)
  if pos.iloc[0,0]=='barcode':pos=pos.iloc[1:]
  pos.columns=['barcode','in_tissue','array_row','array_col','y','x'];pos=pos.set_index('barcode');assert pos.index.is_unique;missing=len(set(bc)-set(pos.index));assert missing==0
  pos=pos.loc[bc].apply(pd.to_numeric);scalep=find('scalefactors');sf=json.loads(gzip.decompress(scalep.read_bytes()) if scalep.suffix=='.gz' else scalep.read_bytes());ip=find('tissue_hires_image')
  if ip.suffix=='.gz':
   import io as bio
   im=Image.open(bio.BytesIO(gzip.decompress(ip.read_bytes()))).convert('RGB')
  else:im=Image.open(ip).convert('RGB')
  im.save(O/(gsm+'_HE.png'));w,h=im.size;xy=pos[['x','y']].to_numpy()*sf['tissue_hires_scalef'];inside=(xy[:,0]>=0)&(xy[:,0]<w)&(xy[:,1]>=0)&(xy[:,1]<h);tissue=pos.in_tissue.to_numpy()==1
  genes={}
  for gene in ['CGNL1','PCSK5','PECAM1','VWF','KDR','EMCN','CDH5','ALB','EPCAM']:
   ix=np.flatnonzero(sy==gene);val=np.asarray(X[ix].sum(axis=0)).ravel() if len(ix) else np.zeros(len(bc));genes[gene]={'feature_present':bool(len(ix)),'positive_tissue_spots':int(((val>0)&tissue).sum()),'total_counts_tissue':float(val[tissue].sum())}
  cgix=np.flatnonzero(sy=='CGNL1');cg=np.asarray(X[cgix].sum(axis=0)).ravel() if len(cgix) else np.zeros(len(bc));z=pos.copy();z['CGNL1_counts']=cg;z['hires_x']=xy[:,0];z['hires_y']=xy[:,1];z.to_csv(O/(gsm+'_spot_linkage.tsv.gz'),sep='\t')
  row={'dataset':g,'sample':gsm,'title':s['title'],'matrix_genes':len(sy),'matrix_spots':len(bc),'tissue_spots':int(tissue.sum()),'CGNL1_feature_present':genes['CGNL1']['feature_present'],'CGNL1_positive_tissue_spots':genes['CGNL1']['positive_tissue_spots'],'CGNL1_total_counts':genes['CGNL1']['total_counts_tissue'],'HE_width':w,'HE_height':h,'coordinate_match_fraction':1-missing/len(bc),'tissue_coordinates_inside_image':float(inside[tissue].mean()),'gene_coverage':genes,'characteristics':s.get('characteristics',[])};rows.append(row)
  ax=axes.ravel()[k];ax.imshow(im);ax.scatter(xy[tissue,0],xy[tissue,1],s=.5,c='cyan',alpha=.15);ax.set_title(g+' / '+gsm+'\nCGNL1 detected: '+str(row['CGNL1_positive_tissue_spots'])+'/'+str(row['tissue_spots']),fontsize=9);ax.axis('off');k+=1
  print(gsm,row['matrix_spots'],row['CGNL1_positive_tissue_spots'],row['tissue_coordinates_inside_image'],flush=True)
(A/'SPATIAL_FILE_AUDIT.json').write_text(json.dumps(rows,indent=2));pd.DataFrame([{k:v for k,v in r.items() if k not in ['gene_coverage','characteristics']} for r in rows]).to_csv(O/'spatial_file_audit.tsv',sep='\t',index=False);fig.tight_layout();fig.savefig(O/'HE_coordinate_contact_sheet.png',dpi=150);plt.close(fig);(A/'SPATIAL_AUDIT_EXIT.json').write_text(json.dumps({'exit_code':0,'sections':len(rows)}))
