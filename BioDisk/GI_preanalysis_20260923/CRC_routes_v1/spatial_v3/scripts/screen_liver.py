import pathlib,json,tarfile,gzip,sys,io,math
import numpy as np,pandas as pd
from scipy import io as sio,stats
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
O=pathlib.Path(sys.argv[1]);B=O.parent;R=O/('results/screen_pilot' if len(sys.argv)>2 else 'results/screen');R.mkdir(exist_ok=True);(R/'figures').mkdir(exist_ok=True)
plan=json.loads((O/'SCREEN_PLAN.json').read_text());genes=plan['candidates'];programs=plan['programs'];jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text());annpath=B/'results/spatial_author_metadata.tsv.gz'
if not annpath.exists():annpath=B/'metadata/spatial_author_metadata.tsv.gz'
jobs=[j for j in jobs if len(sys.argv)<=2 or j['sample']==sys.argv[2]]
ann=pd.read_csv(annpath,sep='\t');ann['sample']=ann.Barcode.str.split('~').str[0];ann['barcode']=ann.Barcode.str.split('~').str[1]
effects=[];regions=[];qcs=[];elig=[];allspots=[];pcoverage=[]
for job in jobs:
 name=job['sample'];patient=job['patient'];rec='Recurrence' in job['sample_type'];path=O/'raw/liver'/(name+'.tar.gz')
 with tarfile.open(path) as tar:
  names=tar.getnames()
  def member(suffix):
   hit=[x for x in names if x.endswith(suffix)];assert len(hit)==1,(name,suffix,hit);return hit[0]
  def readgz(suffix):return gzip.GzipFile(fileobj=tar.extractfile(member(suffix)))
  ft=pd.read_csv(readgz('filtered_feature_bc_matrix/features.tsv.gz'),sep='\t',header=None);sy=ft[1].astype(str).to_numpy();ids=ft[0].astype(str).to_numpy();bc=pd.read_csv(readgz('filtered_feature_bc_matrix/barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str).to_numpy()
  with threadpool_limits(limits=2):X=sio.mmread(readgz('filtered_feature_bc_matrix/matrix.mtx.gz')).tocsr()
  X.eliminate_zeros();assert X.shape==(len(sy),len(bc));assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data));assert pd.Index(ids).is_unique
  a=ann[ann['sample']==name].copy().reset_index(drop=True);assert a.barcode.is_unique;ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all();X=X[:,ix];lib=np.asarray(X.sum(axis=0)).ravel();nfeat=np.asarray(X.getnnz(axis=0)).ravel();bar=a.barcode.to_numpy().copy()
  posmembers=[x for x in names if x.endswith('tissue_positions.csv')];posmembers=posmembers or [x for x in names if x.endswith('tissue_positions_list.csv')];assert len(posmembers)==1
  pos=pd.read_csv(tar.extractfile(posmembers[0]),header=0 if posmembers[0].endswith('tissue_positions.csv') else None)
  if 'barcode' not in pos.columns:pos.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
  a=a.merge(pos,on='barcode',sort=False,validate='one_to_one');assert np.array_equal(a.barcode.to_numpy(),bar)
  sf=json.load(tar.extractfile(member('spatial/scalefactors_json.json')));factor=55/float(sf['spot_diameter_fullres']);plotxy=a[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();image_pitch=np.median(cKDTree(plotxy*factor).query(plotxy*factor,k=2)[0][:,1]);xy=np.c_[50*a.array_col.to_numpy(),50*np.sqrt(3)*a.array_row.to_numpy()];assert ((a.array_col%2)==(a.array_row%2)).all();pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);assert 95<pitch<105,(name,pitch)
  good=(lib>=plan['QC']['UMI_min'])&(nfeat>=plan['QC']['genes_min']);a['UMI']=lib;a['n_genes']=nfeat;a['QC_pass']=good;a['patient']=patient;a['recurrence']=rec;a['x_um']=xy[:,0];a['y_um']=xy[:,1]
  hepatic=(a.Layer3=='Hepatic lobule').to_numpy();macro=(a.Layer3=='Liver macrometastasis tumor').to_numpy();stroma=(a.Layer3=='Liver macrometastasis stroma').to_numpy();micro=(a.Layer3=='Liver micrometastasis tumor').to_numpy()
  a['distance_to_hepatic_um']=cKDTree(xy[hepatic]).query(xy)[0] if hepatic.any() else np.nan
  a['distance_to_macro_um']=cKDTree(xy[macro]).query(xy)[0] if macro.any() else np.nan
  gidx={}
  for g in genes:
   jj=np.where(sy==g)[0];assert len(jj)<=1
   if len(jj):gidx[g]=jj[0];v=X[jj[0]].toarray().ravel();a[g+'_counts']=v;a[g+'_log2CPM']=np.log2(1+v/np.maximum(lib,1)*1e6)
  for pn,gs in programs.items():
   jj=np.where(np.isin(sy,gs))[0];pcoverage.append({'sample':name,'program':pn,'covered':len(jj),'expected':len(gs),'measured':';'.join(sy[jj]),'missing':';'.join(set(gs)-set(sy[jj]))});assert len(jj)>=.7*len(gs),(name,pn,sy[jj]);a[pn+'_score']=np.log2(1+X[jj].toarray()/np.maximum(lib,1)*1e6).mean(axis=0)
  def add_contrast(label,m1,m0,kind):
   n1=int(m1.sum());n0=int(m0.sum());ok=n1>=20 and n0>=20;elig.append({'sample':name,'patient':patient,'recurrence':rec,'contrast':label,'group1_spots':n1,'group0_spots':n0,'eligible':ok})
   if not ok:return
   l1=lib[m1].sum();l0=lib[m0].sum();use=m1|m0
   if kind=='distance':
    z=np.column_stack([a.loc[use,p+'_score'].to_numpy() for p in ['hepatic','stromal','epithelial','endothelial']]+[np.log1p(lib[use])]);std=z.std(axis=0);z=(z-z.mean(axis=0))/np.where(std>0,std,1);design=np.column_stack([np.ones(use.sum()),m1[use].astype(float),z]);rank=np.linalg.matrix_rank(design);condition=np.linalg.cond(design);nuis=np.column_stack([np.ones(use.sum()),z]);pred=nuis@np.linalg.lstsq(nuis,m1[use].astype(float),rcond=None)[0];r2=1-np.sum((m1[use]-pred)**2)/np.sum((m1[use]-np.mean(m1[use]))**2);vif=1/max(1-r2,1e-12)
   for g,i in gidx.items():
    v=X[i].toarray().ravel();y=a[g+'_log2CPM'].to_numpy();common={'sample':name,'patient':patient,'recurrence':rec,'contrast':label,'kind':kind,'gene':g,'n1':n1,'n0':n0,'group1_CPM':v[m1].sum()/l1*1e6,'group0_CPM':v[m0].sum()/l0*1e6,'group1_detected_pct':100*np.mean(v[m1]>0),'group0_detected_pct':100*np.mean(v[m0]>0),'hepatic_score_delta':float(a.loc[m1,'hepatic_score'].mean()-a.loc[m0,'hepatic_score'].mean()),'stromal_score_delta':float(a.loc[m1,'stromal_score'].mean()-a.loc[m0,'stromal_score'].mean()),'epithelial_score_delta':float(a.loc[m1,'epithelial_score'].mean()-a.loc[m0,'epithelial_score'].mean())}
    delta=np.log2(1+common['group1_CPM'])-np.log2(1+common['group0_CPM']);effects.append(dict(common,method='pseudobulk',effect=float(delta)))
    if kind=='distance':
     effects.append(dict(common,method='spot_mean',effect=float(y[m1].mean()-y[m0].mean())));coef=np.linalg.lstsq(design,y[use],rcond=None)[0];valid=rank==design.shape[1] and condition<1e8 and vif<20;effects.append(dict(common,method='score_adjusted',effect=float(coef[1]) if valid else np.nan,near_VIF=vif,condition_number=condition,regression_valid=valid))
  for near,far in [[300,600]]+plan['distance_sensitivities']:
   d=a.distance_to_hepatic_um.to_numpy();add_contrast(f'tumor_near{near}_vs_far{far}',macro&good&(d<=near+1e-6),macro&good&(d>=far-1e-6),'distance')
  for label,m1,m0 in [('hepatic_lobule_vs_macro_tumor',hepatic,macro),('macro_stroma_vs_macro_tumor',stroma,macro),('micro_tumor_vs_macro_tumor',micro,macro)]:add_contrast(label,m1&good,m0&good,'source')
  for lab,aa in a[good].groupby('Layer3'):
   k=aa.index.to_numpy();l=lib[k].sum()
   for g,i in gidx.items():
    v=X[i,k].toarray().ravel();regions.append({'sample':name,'patient':patient,'recurrence':rec,'region':lab,'gene':g,'spots':len(k),'gene_counts':int(v.sum()),'total_UMI':int(l),'CPM':float(v.sum()/l*1e6),'log2CPM':float(np.log2(1+v.sum()/l*1e6)),'detected_pct':float(100*np.mean(v>0)),'mean_spot_log2CPM':float(a.loc[k,g+'_log2CPM'].mean())})
  q={'sample':name,'patient':patient,'recurrence':rec,'matrix_spots':len(bc),'author_spots':len(a),'QC_spots':int(good.sum()),'macro_QC':int((macro&good).sum()),'hepatic_QC':int((hepatic&good).sum()),'micro_QC':int((micro&good).sum()),'median_pitch_um':float(pitch),'image_based_pitch_um_not_used':float(image_pitch),'CGNL1_counts_QC':int(a.loc[good,'CGNL1_counts'].sum()),'CGNL1_detected_pct_QC':float(100*(a.loc[good,'CGNL1_counts']>0).mean())};qcs.append(q)
  # Save compact target/proxy/annotation matrix for targeted reanalysis without large raw loads.
  allspots.append(a)
  fig,axs=plt.subplots(1,3,figsize=(14,4.8));pal={'Hepatic lobule':'#6BAED6','Liver macrometastasis tumor':'#C44E52','Liver macrometastasis stroma':'#999999','Liver micrometastasis tumor':'#AD3BA4','Liver micrometastasis stroma':'#CCA2CD'}
  hires=[x for x in names if x.endswith('spatial/tissue_hires_image.png')]
  if hires:
   im=Image.open(tar.extractfile(hires[0]));axs[0].imshow(im);s=sf['tissue_hires_scalef'];axs[0].scatter(a.pxl_col_in_fullres*s,a.pxl_row_in_fullres*s,c=[pal.get(z,'#DBCE9A') for z in a.Layer3],s=1,alpha=.25)
  else:axs[0].scatter(plotxy[:,0],plotxy[:,1],c=[pal.get(z,'#DBCE9A') for z in a.Layer3],s=2);axs[0].invert_yaxis()
  axs[0].set_title('Pathology: tumor red / hepatic blue / stroma grey',fontsize=9)
  sc=axs[1].scatter(plotxy[:,0],plotxy[:,1],c=np.where(good,a.CGNL1_log2CPM,np.nan),s=2,cmap='viridis',vmin=0,vmax=12);axs[1].invert_yaxis();axs[1].set_title('CGNL1 log2(CPM+1), QC spots',fontsize=10);fig.colorbar(sc,ax=axs[1],shrink=.65)
  d=a.distance_to_hepatic_um.to_numpy();color=np.full(len(a),'#E5E5E5',dtype=object);color[hepatic]='#A6CEE3';color[macro&good&(d<=300+1e-6)]='#E66101';color[macro&good&(d>=600-1e-6)]='#5E3C99';axs[2].scatter(plotxy[:,0],plotxy[:,1],c=color,s=2);axs[2].invert_yaxis();axs[2].set_title('Tumor: near <=300um orange / far >=600um purple',fontsize=9)
  for ax in axs:ax.set_aspect('equal');ax.axis('off')
  fig.suptitle(f'Patient {patient} | {name}'+(' | recurrence' if rec else ' | non-recurrence'),fontsize=12);fig.tight_layout();fig.savefig(R/'figures'/(name+'.png'),dpi=140);plt.close(fig)
  print(json.dumps(q),flush=True)
  # Incremental checkpoint after each completed section.
  pd.DataFrame(qcs).to_csv(R/'sample_QC.tsv',sep='\t',index=False)
pd.DataFrame(regions).to_csv(R/'region_expression.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(R/'section_effects.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(R/'contrast_eligibility.tsv',sep='\t',index=False);pd.DataFrame(pcoverage).to_csv(R/'program_coverage.tsv',sep='\t',index=False);pd.concat(allspots,ignore_index=True).to_csv(R/'spot_target_coordinates.tsv.gz',sep='\t',index=False)
e=pd.DataFrame(effects);patient_rows=[]
for scope in ['baseline','all_timepoints']:
 z=e[~e.recurrence] if scope=='baseline' else e
 for key,v in z.groupby(['patient','gene','contrast','method']):
  val=v.effect.dropna()
  if len(val):patient_rows.append(dict(zip(['patient','gene','contrast','method'],key),scope=scope,effect=float(val.mean()),sections=int(len(val))))
p=pd.DataFrame(patient_rows);p.to_csv(R/'patient_effects.tsv',sep='\t',index=False)
def family(scope,contrast,method,gene):
 if method!='pseudobulk':return scope+'_score_sensitivity'
 if contrast.startswith('tumor_'):
  if scope=='baseline' and contrast=='tumor_near300_vs_far600':return 'single_primary' if gene=='CGNL1' else 'secondary_genes3'
  return scope+'_distance_sensitivity'
 return scope+'_source_comparisons'
rows=[]
for key,z in p.groupby(['scope','gene','contrast','method']):
 vals=z.effect.to_numpy();n=len(vals);d=dict(zip(['scope','gene','contrast','method'],key));d.update(n_patients=n,positive=int((vals>0).sum()),negative=int((vals<0).sum()),mean_effect=float(vals.mean()),median_effect=float(np.median(vals)),LOO_min=float(min(np.delete(vals,i).mean() for i in range(n))) if n>1 else np.nan,LOO_max=float(max(np.delete(vals,i).mean() for i in range(n))) if n>1 else np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  se=stats.sem(vals);crit=stats.t.ppf(.975,n-1);d.update(P=float(stats.ttest_1samp(vals,0).pvalue),CI_low=float(vals.mean()-crit*se),CI_high=float(vals.mean()+crit*se))
 d['family']=family(d['scope'],d['contrast'],d['method'],d['gene']);rows.append(d)
s=pd.DataFrame(rows);s['FDR']=np.nan
expected={'secondary_genes3':3,'baseline_distance_sensitivity':8,'all_timepoints_distance_sensitivity':12,'baseline_source_comparisons':12,'all_timepoints_source_comparisons':12,'baseline_score_sensitivity':24,'all_timepoints_score_sensitivity':24}
for fam,zz in s.groupby('family'):
 if fam=='single_primary':continue
 N=expected[fam];assert len(zz)<=N;pv=np.r_[zz.P.fillna(1).to_numpy(),np.ones(N-len(zz))];order=np.argsort(pv);adj=np.minimum.accumulate((pv[order]*N/np.arange(1,N+1))[::-1])[::-1];q=np.empty(N);q[order]=np.minimum(adj,1);q=q[:len(zz)];q[zz.P.isna().to_numpy()]=np.nan;s.loc[zz.index,'FDR']=q
s.to_csv(R/'summary.tsv',sep='\t',index=False);(R/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'sections':len(qcs),'patients':len(set(x['patient'] for x in qcs)),'baseline_sections':sum(not x['recurrence'] for x in qcs),'QC_spots':sum(x['QC_spots'] for x in qcs),'primary':s[(s.family=='single_primary')].to_dict('records'),'scope':'Pathology-defined mixed spots; no cell-specific causal inference; patient-level effects only.'},indent=2));print(s[(s.scope=='baseline')&(s.gene=='CGNL1')].to_string(index=False));print('SPATIAL_SCREEN_COMPLETE',flush=True)
