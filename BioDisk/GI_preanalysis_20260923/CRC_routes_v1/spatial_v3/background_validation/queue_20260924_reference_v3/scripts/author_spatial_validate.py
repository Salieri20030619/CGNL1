from common import *
from scipy.spatial import cKDTree
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Q/'V07_AUTHOR_SPATIAL';D.mkdir(exist_ok=True);A=Q/'AUTHOR_VALIDATION'
if (A/'author_counts.npz').exists():X=sparse.load_npz(A/'author_counts.npz')
else:
 with threadpool_limits(limits=2):X=io.mmread(A/'author_counts.mtx').tocsr()
sy=pd.read_csv(A/'author_genes.tsv',sep='\t').gene.astype(str).to_numpy();v=pd.read_csv(A/'author_spots.tsv',sep='\t');bar=pd.read_csv(A/'author_barcodes.tsv',sep='\t');assert np.array_equal(v.batch,bar.batch);v['barcode']=bar.barcode;assert X.shape==(len(sy),len(v));assert (X.data>=0).all() and (X.data==np.floor(X.data)).all()
lib=np.asarray(X.sum(axis=0)).ravel();assert np.array_equal(lib,v.library_provided_genes)
sparse.save_npz(A/'author_counts.npz',X)
qc=[];effects=[];checks=[]
for patient,t in v.groupby('batch',sort=True):
 ix=t.index.to_numpy();Y=X[:,ix];L=lib[ix];xy=t[['imagecol','imagerow']].to_numpy();rawdist=cKDTree(xy).query(xy,k=2)[0][:,1];pitch=np.median(rawdist);cv=np.std(rawdist)/pitch
 if cv>=.1:
  qc.append({'patient':int(patient),'all_spots':len(t),'status':'coordinate_calibration_not_eligible','raw_coordinate_pitch':float(pitch),'nearest_pitch_CV':float(cv)});continue
 xy=xy/pitch*100
 if patient in [4,5]:
  name='patient'+str(patient);z=np.load(Q/'V04_INPUT_QC'/(name+'_features.npz'));Z=sparse.load_npz(Q/'V04_INPUT_QC'/(name+'_counts.npz'));bi=pd.Index(z['barcodes']).get_indexer(t.barcode);rawgenes=pd.Index(z['genes']);uniq=~rawgenes.duplicated(keep=False);index={g:j for j,g in enumerate(rawgenes) if uniq[j]};present=np.array([g in index for g in sy]);gi=np.array([index[g] for g in sy[present]]);assert (bi>=0).all() and present.sum()>20000;assert all(g in index for g in C['targets']);diff=Y[present]-Z[gi][:,bi];diff.eliminate_zeros();from scipy.stats import spearmanr;rho=float(spearmanr(np.asarray(Y[present].sum(axis=0)).ravel(),np.asarray(Z[gi][:,bi].sum(axis=0)).ravel()).statistic);assert rho>.99,'sample alignment not supported';checks.append({'patient':int(patient),'matched_spots':len(t),'verified_unique_genes':int(present.sum()),'excluded_ambiguous_or_absent_symbols':int((~present).sum()),'all_verified_counts_identical':bool(diff.nnz==0),'mismatched_entries':int(diff.nnz),'shared_library_spearman':rho,'input_classification':'author processed integer counts; not original raw UMI','median_retained_UMI_fraction':float(np.median(L/np.asarray(Z[:,bi].sum(axis=0)).ravel()))})

 nct=t[['hepatocytes+cholangiocytes','tumor_epithelial','myeloid','stromal']].to_numpy();frac=nct/np.maximum(nct.sum(axis=1,keepdims=True),1e-12);good=(L>=500)&(Y.getnnz(axis=0)>=200);hephits=Y[np.isin(sy,['ALB','CPS1','APOA1','APOA2','ASGR1','TTR','HAMP'])].getnnz(axis=0);hep=good&(frac[:,0]>=.7)&(hephits>=2);tum=good&(frac[:,1]>=.6)&(frac[:,0]<.1)&(frac[:,3]<.2)
 rec={'patient':int(patient),'all_spots':len(t),'QC_spots':int(good.sum()),'hepatic_spots':int(hep.sum()),'tumor_like_spots':int(tum.sum()),'raw_coordinate_pitch':float(pitch),'nearest_pitch_CV':float(cv)}
 if hep.sum()<20 or tum.sum()<20:rec['status']='insufficient_regions';qc.append(rec);continue
 dist=cKDTree(xy[hep]).query(xy)[0];dd,jj=cKDTree(xy).query(xy,k=7,distance_upper_bound=110);valid=np.isfinite(dd[:,1:])&(jj[:,1:]<len(t));neighbors=tum[np.minimum(jj[:,1:],len(t)-1)]&valid;strict=tum&(valid.sum(axis=1)==6)&(neighbors.sum(axis=1)>=5);groups={'strict5':(strict&(dist>110)&(dist<=300),strict&(dist>=600)),'original_region':(tum&(dist<=300),tum&(dist>=600))}
 rec.update(status='regions_available',strict_near=int(groups['strict5'][0].sum()),strict_far=int(groups['strict5'][1].sum()))
 for mask,(near,far) in groups.items():
  for gene in C['targets']:
   ids=np.flatnonzero(sy==gene)
   if len(ids)!=1 or near.sum()<20 or far.sum()<20:continue
   counts=Y[ids[0]].toarray().ravel();nearCPM=counts[near].sum()/L[near].sum()*1e6;farCPM=counts[far].sum()/L[far].sum()*1e6;effects.append({'patient':int(patient),'gene':gene,'mask':mask,'effect':float(np.log2(1+nearCPM)-np.log2(1+farCPM)),'near_spots':int(near.sum()),'far_spots':int(far.sum()),'near_CPM':float(nearCPM),'far_CPM':float(farCPM)})
 q=t[['barcode','imagerow','imagecol']].copy();q['distance_um']=dist;q['hepatic']=hep;q['tumor_like']=tum;q['strict_near']=groups['strict5'][0];q['strict_far']=groups['strict5'][1]
 for gene in C['targets']:
  g=np.flatnonzero(sy==gene)
  if len(g)==1:q[gene+'_counts']=Y[g[0]].toarray().ravel()
 q.to_csv(D/('patient'+str(patient)+'_regions.tsv.gz'),sep='\t',index=False)
 fig,axs=plt.subplots(1,3,figsize=(12,4));color=np.full(len(t),'#DDD',dtype=object);color[hep]='#619EC4';color[groups['strict5'][0]]='#C46D43';color[groups['strict5'][1]]='#7963A8';axs[0].scatter(xy[:,0],xy[:,1],c=color,s=4);axs[0].set_title('Hepatic / near / far regions')
 for ax,g in zip(axs[1:],['CGNL1','SLC9B2']):
  i=np.flatnonzero(sy==g)
  if len(i)==1:im=ax.scatter(xy[:,0],xy[:,1],c=np.log2(1+Y[i[0]].toarray().ravel()/np.maximum(L,1)*1e6),s=4,cmap='viridis');fig.colorbar(im,ax=ax,shrink=.7)
  ax.set_title(g+' log2(CPM+1)')
 for ax in axs:ax.invert_yaxis();ax.set_aspect('equal');ax.axis('off')
 fig.suptitle('Author batch '+str(patient)+' | descriptive validation');fig.tight_layout();fig.savefig(D/('patient'+str(patient)+'.png'),dpi=150);plt.close(fig);qc.append(rec)
pd.DataFrame(qc).to_csv(D/'QC.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(D/'patient_effects.tsv',sep='\t',index=False);dump(D/'COUNT_INTEGRITY.json',checks)
rows=[]
if effects:
 e=pd.DataFrame(effects)
 for (g,m),t in e.groupby(['gene','mask']):rows.append({'gene':g,'mask':m,'patients':len(t),'positive':int((t.effect>0).sum()),'mean_effect':float(t.effect.mean()),'median_effect':float(t.effect.median())})
dump(D/'RESULT.json',{'status':'exploratory_external_computed','QC':qc,'effects':rows,'source':'author integrated NicheDE object with RCTD composition; barcodes/genes mapped; small count discrepancies explicitly retained, no raw-UMI equality claim','scope':'5 external CRC liver metastasis patients, n<10 descriptive only. Author processed integer counts: no inference as original raw UMI; findings require raw-data confirmation. 21708 provided genes for CPM denominator; bulk scales not pooled. Author hepatic group includes cholangiocytes; >=2 hepatocyte marker detections required. Pixel spacing calibrated per sample; 110um exclusion handles pitch rounding. Author regions not independent malignancy confirmation.'})
