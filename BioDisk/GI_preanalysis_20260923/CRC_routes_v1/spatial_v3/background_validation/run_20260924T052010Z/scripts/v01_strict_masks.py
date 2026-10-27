import pathlib,sys,json,hashlib
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=pathlib.Path(sys.argv[1]);C=json.loads((R/'CONFIG.json').read_text());S=pathlib.Path(C['remote_spatial_root']);D=R/'V01';D.mkdir(exist_ok=True);(D/'figures').mkdir(exist_ok=True)
input_path=S/'results/screen/spot_target_coordinates.tsv.gz';a=pd.read_csv(input_path,sep='\t');assert not a.duplicated(['sample','barcode']).any();a=a[~a.recurrence].copy();effects=[];elig=[];masks=[]
for name,v in a.groupby('sample',sort=True):
 v=v.reset_index(drop=True);xy=v[['x_um','y_um']].to_numpy();dist,ix=cKDTree(xy).query(xy,k=7,distance_upper_bound=100+1e-4);dist=dist[:,1:];ix=ix[:,1:];valid=np.isfinite(dist)&(ix<len(v));macro=(v.Layer3=='Liver macrometastasis tumor').to_numpy();neighbours=macro[np.minimum(ix,len(v)-1)]&valid;nt=neighbours.sum(axis=1);nv=valid.sum(axis=1);good=v.QC_pass.to_numpy();dh=v.distance_to_hepatic_um.to_numpy();fig,axs=plt.subplots(1,2,figsize=(10,4.5));patient=int(v.patient.iloc[0])
 for purity,ax in zip([5,6],axs):
  keep=macro&good&(nv==6)&(nt>=purity);near=keep&(dh>100+1e-6)&(dh<=300+1e-6);far=keep&(dh>=600-1e-6);n1=int(near.sum());n0=int(far.sum());rule=str(purity)+'_of_6';ok=n1>=20 and n0>=20;elig.append({'sample':name,'patient':patient,'rule':rule,'near_spots':n1,'far_spots':n0,'eligible':ok})
  masks.append(pd.DataFrame({'sample':name,'barcode':v.barcode,'patient':patient,'rule':rule,'six_neighbours_annotated':nv==6,'tumour_neighbours':nt,'strict_tumour':keep,'near':near,'far':far}))
  cc=np.full(len(v),'#DEDEDE',dtype=object);cc[near]='#D06B45';cc[far]='#6750A1';px=v[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();ax.scatter(px[:,0],px[:,1],c=cc,s=2);ax.invert_yaxis();ax.set_aspect('equal');ax.axis('off');ax.set_title(f'{rule}: near={n1}, far={n0}',fontsize=10)
  if not ok:continue
  use=near|far;lib=v.UMI.to_numpy();z=np.column_stack([v.loc[use,p].to_numpy() for p in ['hepatic_score','stromal_score','epithelial_score','endothelial_score']]+[np.log1p(lib[use])]);sd=z.std(axis=0);z=(z-z.mean(axis=0))/np.where(sd>0,sd,1);X=np.c_[np.ones(use.sum()),near[use].astype(float),z];N=np.c_[np.ones(use.sum()),z];pred=N@np.linalg.lstsq(N,near[use].astype(float),rcond=None)[0];r2=1-np.sum((near[use]-pred)**2)/np.sum((near[use]-near[use].mean())**2);vif=1/max(1-r2,1e-12);regok=np.linalg.matrix_rank(X)==X.shape[1] and np.linalg.cond(X)<1e8 and vif<20
  for gene in C['genes']:
   counts=v[gene+'_counts'].to_numpy();y=np.log2(1+counts/np.maximum(lib,1)*1e6);pb=np.log2(1+counts[near].sum()/lib[near].sum()*1e6)-np.log2(1+counts[far].sum()/lib[far].sum()*1e6);base={'sample':name,'patient':patient,'rule':rule,'gene':gene,'near_spots':n1,'far_spots':n0,'near_VIF':vif}
   for method,value in [('pseudobulk',pb),('spot_mean',y[near].mean()-y[far].mean()),('score_adjusted',np.linalg.lstsq(X,y[use],rcond=None)[0][1] if regok else np.nan)]:effects.append(dict(base,method=method,effect=float(value)))
 fig.suptitle(f'Patient {patient} | {name} | Internal robustness');fig.tight_layout();fig.savefig(D/'figures'/(name+'.png'),dpi=130);plt.close(fig);print(name,'processed',flush=True)
pd.DataFrame(elig).to_csv(D/'eligibility.tsv',sep='\t',index=False);pd.concat(masks,ignore_index=True).to_csv(D/'masks.tsv.gz',sep='\t',index=False)
if not effects:(D/'RESULT_STATUS.json').write_text(json.dumps({'status':'insufficient_data','reason':'No section met20 spots per group'}));sys.exit(0)
e=pd.DataFrame(effects);e.to_csv(D/'section_effects.tsv',sep='\t',index=False);p=e.groupby(['patient','rule','gene','method']).effect.mean().dropna().reset_index();p.to_csv(D/'patient_effects.tsv',sep='\t',index=False);rows=[]
for key,v in p.groupby(['rule','gene','method']):
 x=v.effect.to_numpy();n=len(x);r=dict(zip(['rule','gene','method'],key));r.update(n_patients=n,positive=int((x>0).sum()),negative=int((x<0).sum()),mean_effect=float(x.mean()),median_effect=float(np.median(x)),P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  se=stats.sem(x);t=stats.t.ppf(.975,n-1);r.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=float(x.mean()-t*se),CI_high=float(x.mean()+t*se))
 r['family']='score_sensitivity' if key[2]!='pseudobulk' else ('strict6_sensitivity' if key[0]=='6_of_6' else ('single_primary' if key[1]=='CGNL1' else 'secondary_genes'));rows.append(r)
s=pd.DataFrame(rows);s['FDR']=np.nan
for fam,z in s.groupby('family'):
 if fam=='single_primary':continue
 n=C['V01']['families'][fam];assert len(z)<=n;v=np.r_[z.P.fillna(1),np.ones(n-len(z))];o=np.argsort(v);q=np.minimum.accumulate((v[o]*n/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[o]=np.minimum(q,1);out=out[:len(z)];out[z.P.isna().to_numpy()]=np.nan;s.loc[z.index,'FDR']=out
s.to_csv(D/'summary.tsv',sep='\t',index=False)
with input_path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
(D/'RESULT_STATUS.json').write_text(json.dumps({'status':'computed_pending_review','input_sha256':sha,'baseline_sections':a['sample'].nunique(),'primary':s[s.family=='single_primary'].to_dict('records'),'scope':'Internal spatial robustness; mixed spots remain possible; not independent replication or causal evidence.'},indent=2));print('V01 outputs ready',flush=True)
