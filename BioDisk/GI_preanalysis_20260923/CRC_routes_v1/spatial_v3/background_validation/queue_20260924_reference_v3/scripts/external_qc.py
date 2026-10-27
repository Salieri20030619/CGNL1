from common import *
import zipfile,io as bio,h5py
D=Q/'V04_INPUT_QC';D.mkdir(exist_ok=True);rows=[]
for sample,archive,root in [('patient4','22261_visium.zip','22261_visium/'),('patient5','22996_visium.zip','22996_visium/outs/')]:
 with zipfile.ZipFile(I/'inputs/NicheDE'/archive) as z:
  names=z.namelist();hpath=root+'filtered_feature_bc_matrix.h5';assert hpath in names
  with h5py.File(bio.BytesIO(z.read(hpath)),'r') as h:
   m=h['matrix'];X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]).tocsr();sy=np.array([x.decode() for x in m['features/name'][:]]);bc=np.array([x.decode() for x in m['barcodes'][:]])
  pos=pd.read_csv(bio.BytesIO(z.read(root+'spatial/tissue_positions_list.csv')),header=None);pos.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'];pos=pos.set_index('barcode').reindex(bc);assert pos.in_tissue.notna().all();pos.index.name='barcode';pos.to_csv(D/(sample+'_positions.tsv'),sep='\t');assert X.shape==(len(sy),len(bc));sparse.save_npz(D/(sample+'_counts.npz'),X);np.savez_compressed(D/(sample+'_features.npz'),genes=sy.astype(str),barcodes=bc.astype(str));(D/(sample+'_scalefactors.json')).write_bytes(z.read(root+'spatial/scalefactors_json.json'))
  for im in ['tissue_hires_image.png','tissue_lowres_image.png']:
   p=root+'spatial/'+im
   if p in names:(D/(sample+'_'+im)).write_bytes(z.read(p))
  quarantined=[p for p in names if 'Breast_Cancer' in p and not p.startswith('__MACOSX')];lib=np.asarray(X.sum(axis=0)).ravel();good=(lib>=500)&(np.asarray(X.getnnz(axis=0)).ravel()>=200);g=np.flatnonzero(sy=='CGNL1');rows.append({'sample':sample,'matrix_root':root,'spots':len(bc),'QC_spots':int(good.sum()),'CGNL1_features':len(g),'CGNL1_QC_detected_pct':float((X[g][:,good].toarray().ravel()>0).mean()*100) if len(g)==1 else np.nan,'unrelated_example_entries_excluded':len(quarantined),'clinical_site':'CRC; exact colon/rectum not confirmed','morphology_review':'pending'});dump(D/(sample+'_excluded_entries.json'),quarantined)
pd.DataFrame(rows).to_csv(D/'QC.tsv',sep='\t',index=False);dump(D/'RESULT.json',{'status':'whitelist_matrices_ready','samples':2,'external_validation_status':'not_complete; model and morphology review pending'})
