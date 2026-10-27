from common import *
D=Q/'V08_DISTANCE_RANK';D.mkdir(exist_ok=True);A=Q/'AUTHOR_VALIDATION'
dump(D/'PLAN.json',{'status':'secondary_geometry_adapted_exploration','reason':'Frozen <=300um strict near region had zero spots in all3 hepatic-containing external sections; preserve primary non-estimability','endpoint':'Within strict tumor-like region, lowest vs highest distance quartiles. Masks chosen solely by geometry, before target effects. All4 predeclared genes. >=20 spots/group, report actual distances. No P/FDR for3 patients. Does not test direct contact.'})
v=pd.read_csv(A/'author_spots.tsv',sep='\t');sy=pd.read_csv(A/'author_genes.tsv',sep='\t').gene.to_numpy();X=sparse.load_npz(A/'author_counts.npz');rows=[];qc=[]
from scipy.spatial import cKDTree
for patient in [1,2,3]:
 f=Q/'V07_AUTHOR_SPATIAL'/('patient'+str(patient)+'_regions.tsv.gz');t=pd.read_csv(f,sep='\t');vv=v[v.batch==patient];assert len(t)==len(vv);xy=t[['imagecol','imagerow']].to_numpy();pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);xy=xy/pitch*100;tum=t.tumor_like.to_numpy(bool);dd,ii=cKDTree(xy).query(xy,k=7,distance_upper_bound=110);valid=np.isfinite(dd[:,1:])&(ii[:,1:]<len(t));nn=tum[np.minimum(ii[:,1:],len(t)-1)]&valid;strict=tum&(valid.sum(axis=1)==6)&(nn.sum(axis=1)>=5);dist=t.distance_um.to_numpy();d=dist[strict];q25,q75=np.quantile(d,[.25,.75]);near=strict&(dist<=q25);far=strict&(dist>=q75);Y=X[:,vv.index];lib=np.asarray(Y.sum(axis=0)).ravel();rec={'patient':patient,'strict_spots':int(strict.sum()),'near_spots':int(near.sum()),'far_spots':int(far.sum()),'near_min_um':float(dist[near].min()),'near_max_um':float(dist[near].max()),'far_min_um':float(dist[far].min()),'far_max_um':float(dist[far].max())};qc.append(rec)
 assert near.sum()>=20 and far.sum()>=20
 for g in C['targets']:
  ix=np.flatnonzero(sy==g)
  if len(ix)!=1:continue
  a=Y[ix[0]].toarray().ravel();nearCPM=a[near].sum()/lib[near].sum()*1e6;farCPM=a[far].sum()/lib[far].sum()*1e6;rows.append({**rec,'gene':g,'effect':float(np.log2(1+nearCPM)-np.log2(1+farCPM)),'near_CPM':float(nearCPM),'far_CPM':float(farCPM)})
e=pd.DataFrame(rows);e.to_csv(D/'patient_effects.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(D/'geometry.tsv',sep='\t',index=False);summary=[]
for g,t in e.groupby('gene'):summary.append({'gene':g,'patients':len(t),'positive':int((t.effect>0).sum()),'mean':float(t.effect.mean()),'median':float(t.effect.median())})
dump(D/'RESULT.json',{'status':'secondary_exploratory_computed','summary':summary,'geometry':qc,'scope':'Secondary distance-ranked endpoint after primary geometry failed. Author processed counts. Three patients descriptive only; not independent confirmation of <=300um contact endpoint.'})
