from common import *
from scipy import stats
from scipy.spatial import cKDTree
D=Q/'F02_SPATIAL';D.mkdir(exist_ok=True);prior=pathlib.Path(C['prior_run']);rows=[];elig=[]
def calculate(dataset,sample,patient,counts,lib,dist,strict):
 for quant in [.2,.25,1/3]:
  valid=strict&np.isfinite(dist)
  if valid.sum()<60:elig.append({'dataset':dataset,'sample':sample,'q':quant,'status':'too_few_strict_spots'});continue
  lo,hi=np.quantile(dist[valid],[quant,1-quant]);near=valid&(dist<=lo);far=valid&(dist>=hi);ok=near.sum()>=20 and far.sum()>=20 and lo<hi
  elig.append({'dataset':dataset,'sample':sample,'q':quant,'near_spots':int(near.sum()),'far_spots':int(far.sum()),'near_max_um':float(lo),'far_min_um':float(hi),'status':'eligible' if ok else 'insufficient_groups'})
  if not ok:continue
  for g,z in counts.items():
   a=z[near].sum()/lib[near].sum()*1e6;b=z[far].sum()/lib[far].sum()*1e6
   for pseudo in [.1,1,5]:rows.append({'dataset':dataset,'sample':sample,'patient':str(patient),'gene':g,'q':quant,'pseudocount':pseudo,'effect':float(np.log2(a+pseudo)-np.log2(b+pseudo)),'near_CPM':float(a),'far_CPM':float(b),'near_counts':int(z[near].sum()),'far_counts':int(z[far].sum()),'near_max_um':float(lo),'far_min_um':float(hi)})
a=pd.read_csv(S/'results/screen/spot_target_coordinates.tsv.gz',sep='\t');a=a[~a.recurrence];m=pd.read_csv(I/'V01/masks.tsv.gz',sep='\t');m=m[m.rule=='5_of_6']
for name,v in a.groupby('sample'):
 mm=m[m['sample']==name].set_index('barcode').reindex(v.barcode);strict=mm.strict_tumour.fillna(False).to_numpy(bool)&v.QC_pass.to_numpy(bool);calculate('discovery',name,v.patient.iloc[0],{g:v[g+'_counts'].to_numpy() for g in C['final_targets']},v.UMI.to_numpy(),v.distance_to_hepatic_um.to_numpy(),strict)
A=prior/'AUTHOR_VALIDATION';X=sparse.load_npz(A/'author_counts.npz');sy=pd.read_csv(A/'author_genes.tsv',sep='\t').gene.to_numpy();v=pd.read_csv(A/'author_spots.tsv',sep='\t')
for patient in [1,2,3]:
 t=pd.read_csv(prior/'V07_AUTHOR_SPATIAL'/('patient'+str(patient)+'_regions.tsv.gz'),sep='\t');vv=v[v.batch==patient];xy=t[['imagecol','imagerow']].to_numpy();pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);xy=xy/pitch*100;tum=t.tumor_like.to_numpy(bool);dd,ii=cKDTree(xy).query(xy,k=7,distance_upper_bound=110);valid=np.isfinite(dd[:,1:])&(ii[:,1:]<len(t));nn=tum[np.minimum(ii[:,1:],len(t)-1)]&valid;strict=tum&(valid.sum(axis=1)==6)&(nn.sum(axis=1)>=5);Y=X[:,vv.index];lib=np.asarray(Y.sum(axis=0)).ravel();calculate('external_author_processed','patient'+str(patient),patient,{g:Y[np.flatnonzero(sy==g)[0]].toarray().ravel() for g in C['final_targets']},lib,t.distance_um.to_numpy(),strict)
pd.DataFrame(elig).to_csv(D/'geometry_eligibility.tsv',sep='\t',index=False);e=pd.DataFrame(rows);e.to_csv(D/'section_effects.tsv',sep='\t',index=False);p=e.groupby(['dataset','patient','gene','q','pseudocount'],as_index=False).effect.mean();p.to_csv(D/'patient_effects.tsv',sep='\t',index=False);out=[]
for key,v in p.groupby(['dataset','gene','q','pseudocount']):
 x=v.effect.to_numpy();n=len(x);r=dict(zip(['dataset','gene','q','pseudocount'],key));r.update(n_patients=n,positive=int((x>0).sum()),mean=float(x.mean()),median=float(np.median(x)),P=np.nan,FDR=np.nan,CI_low=np.nan,CI_high=np.nan,LOO_min=float(min(np.delete(x,k).mean() for k in range(n))) if n>1 else np.nan,LOO_max=float(max(np.delete(x,k).mean() for k in range(n))) if n>1 else np.nan)
 if n>=10:
  se=stats.sem(x);t=stats.t.ppf(.975,n-1);r.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=float(x.mean()-t*se),CI_high=float(x.mean()+t*se))
 out.append(r)
s=pd.DataFrame(out)
for key,z in s.groupby(['dataset','q','pseudocount']):
 qs=bh(z.P.to_numpy(),2);qs[z.P.isna().to_numpy()]=np.nan;s.loc[z.index,'FDR']=qs
s.to_csv(D/'summary.tsv',sep='\t',index=False);main=s[np.isclose(s.q,.25)&np.isclose(s.pseudocount,1)]
dump(D/'RESULT.json',{'status':'final_geometry_sensitivity_completed','primary_exploratory_summary':main.to_dict('records'),'scope':'Distance-ranked secondary endpoint, not original <=300um validation. All quantiles/pseudocounts retained. External n3 processed-count data remain exploratory.'})
