from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,gzip,tarfile,io,datetime,itertools
import numpy as np,pandas as pd,h5py
from scipy import sparse,io as sio,stats
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'reassessment_20261008_v1';O=R/'results';spec=json.loads((R/'REASSESSMENT_SPEC.json').read_text());panels=json.loads((B/'v3_pilot_20261006/plan_v3.json').read_text())['panels'];ann=pd.read_csv(B/'results/spatial_author_metadata.tsv.gz',sep='\t');ann['sample']=ann.Barcode.str.split('~').str[0];ann['barcode']=ann.Barcode.str.split('~').str[1]
jobs=[dict(x,cohort='GSE294385') for x in json.loads((B/'spatial_v3/metadata/LIVER_DOWNLOAD_PLAN.json').read_text()) if 'Recurrence' not in x['sample_type']]
fdir=B/'expanded_screen_20260925/spatial_data/GSE217414'
jobs += [dict(sample=p.name.split('_')[0],patient=i+1,cohort='GSE217414',path=str(p)) for i,p in enumerate(sorted(fdir.glob('*_filtered_feature_bc_matrix.h5')))]
def load_job(job):
 sample=job['sample']
 if job['cohort']=='GSE294385':
  path=B/'spatial_v3/raw/liver'/f'{sample}.tar.gz';files={}
  with tarfile.open(path,'r|gz') as tar:
   for m in tar:
    if not m.isfile():continue
    if m.name.endswith('filtered_feature_bc_matrix/features.tsv.gz'):
     with gzip.GzipFile(fileobj=tar.extractfile(m)) as f:genes=pd.read_csv(f,sep='\t',header=None)[1].astype(str).to_numpy()
    elif m.name.endswith('filtered_feature_bc_matrix/barcodes.tsv.gz'):
     with gzip.GzipFile(fileobj=tar.extractfile(m)) as f:bc=pd.read_csv(f,sep='\t',header=None)[0].astype(str).to_numpy()
    elif m.name.endswith('filtered_feature_bc_matrix/matrix.mtx.gz'):
     with gzip.GzipFile(fileobj=tar.extractfile(m)) as f,threadpool_limits(limits=4):X=sio.mmread(f).tocsr()
    elif m.name.endswith('spatial/tissue_positions.csv'):files['pos']=tar.extractfile(m).read()
    elif m.name.endswith('spatial/tissue_positions_list.csv'):files['poslist']=tar.extractfile(m).read()
  a=ann[ann['sample'].eq(sample)].copy().reset_index(drop=True);ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all();X=X[:,ix];X.eliminate_zeros();umi=np.asarray(X.sum(0)).ravel();keep=(umi>=500)&(X.getnnz(0)>=200)
  pos=pd.read_csv(io.BytesIO(files['pos'])) if 'pos' in files else pd.read_csv(io.BytesIO(files['poslist']),header=None,names=['barcode','in_tissue','array_row','array_col','py','px']);pos=pos.set_index('barcode').loc[a.barcode];xy=np.c_[pos.array_col.to_numpy()*50,pos.array_row.to_numpy()*50*np.sqrt(3)];region=a.Layer3.astype(str).to_numpy()
 else:
  with h5py.File(job['path']) as h:
   g=h['matrix'];X=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).tocsr();genes=np.array([x.decode() for x in g['features/name'][:]]);bc=np.array([x.decode() for x in g['barcodes'][:]])
  sd=next((fdir/(sample+'_unpacked')).glob('*spatial'));pos=pd.read_csv(sd/'tissue_positions_list.csv',header=None,names=['barcode','tissue','ar','ac','py','px']).set_index('barcode').loc[bc];sf=json.loads((sd/'scalefactors_json.json').read_text());umi=np.asarray(X.sum(0)).ravel();keep=(umi>=500)&pos.tissue.eq(1).to_numpy();xy=pos[['px','py']].to_numpy()*65/sf['spot_diameter_fullres'];region=np.full(len(umi),'unavailable')
 X=X[:,keep];umi=umi[keep];xy=xy[keep];region=region[keep]
 if len(set(genes))!=len(genes):
  unique,inv=np.unique(genes,return_inverse=True);C=sparse.csr_matrix((np.ones(len(genes)),(inv,np.arange(len(genes)))),shape=(len(unique),len(genes)));X=(C@X).tocsr();genes=unique
 return X,genes,umi,xy,region

# Gene selection uses expression and detection only; no neighborhood or target association is inspected here.
statsrows=[];common=None
for job in jobs:
 print('R1_MATCHING_INPUT',job['sample'],flush=True);X,genes,umi,xy,regions=load_job(job);common=set(genes) if common is None else common&set(genes)
 if job['cohort']=='GSE294385':
  mean=np.asarray(X.multiply(1e4/np.maximum(umi,1)).mean(1)).ravel();detect=np.asarray((X>0).mean(1)).ravel();statsrows.append(pd.DataFrame({'gene':genes,'patient':job['patient'],'sample':job['sample'],'mean_CP10k':mean,'detection':detect}))
 del X
table=pd.concat(statsrows);table.to_csv(O/'R1_expression_matching_inputs.tsv.gz',sep='\t',index=False);pmean=table.groupby(['patient','gene'],as_index=False)[['mean_CP10k','detection']].mean();mean=pmean.groupby('gene')[['mean_CP10k','detection']].mean();mean=mean.loc[sorted(common)].copy();target=mean.loc['CGNL1'];excluded={'CGNL1','ENSG00000128849','JACOP','PCING'}|set(sum(panels.values(),[]))|set(sum(spec['R2_lineage_state']['programs'].values(),[]))
pool=mean[~mean.index.isin(excluded)&mean.mean_CP10k.between(target.mean_CP10k/2,target.mean_CP10k*2)&mean.detection.between(target.detection-.15,target.detection+.15)].copy();assert len(pool)>=20,('expression_matched_controls_not_qualified',len(pool))
def transformed(t):
 p=np.clip(t.detection.to_numpy(),1e-5,1-1e-5);return np.c_[np.log(t.mean_CP10k.to_numpy()+1e-6),np.log(p/(1-p))]
allz=transformed(mean);tz=transformed(target.to_frame().T)[0];scale=np.std(allz,axis=0);pool['distance']=np.linalg.norm((transformed(pool)-tz)/scale,axis=1);pool=pool.sort_values(['distance'],kind='stable');chosen=pool.head(50).copy();chosen.index.name='gene';chosen.reset_index().to_csv(O/'R1_selected_controls.tsv',sep='\t',index=False);frozen={'selected_before_association':True,'eligible_controls':len(pool),'chosen_controls':len(chosen),'target_mean_CP10k':float(target.mean_CP10k),'target_detection':float(target.detection),'controls':chosen.index.tolist(),'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selection_uses_target_association':False};(O/'R1_CONTROL_FREEZE.json').write_text(json.dumps(frozen,indent=2));print('R1_CONTROLS_FROZEN',len(chosen),flush=True)
targets=['CGNL1']+chosen.index.tolist();rows=[];qc=[]
def standard(v):
 s=v.std(0);return (v-v.mean(0))/np.where(s>1e-12,s,1),s
for job in jobs:
 print('R1_MODEL_INPUT',job['sample'],flush=True);X,genes,umi,xy,regions=load_job(job);index={g:i for i,g in enumerate(genes)};wanted=set(targets)|set(sum(panels.values(),[]));values={g:np.log1p(X[index[g]].toarray().ravel()/np.maximum(umi,1)*1e4) for g in wanted if g in index};scores={k:np.mean([values[g] for g in gs if g in values],axis=0) if sum(g in values for g in gs)>=2 else np.full(len(umi),np.nan) for k,gs in panels.items()};Y=np.column_stack([values[g] for g in targets]);tree=cKDTree(xy)
 for radius in [200,100,300]:
  pair=tree.query_pairs(radius+1e-6,output_type='ndarray');rr=np.r_[pair[:,0],pair[:,1]];cc=np.r_[pair[:,1],pair[:,0]];G=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(umi),len(umi)));degree=np.asarray(G.sum(1)).ravel();W=sparse.diags(1/np.maximum(degree,1))@G;neighbor=np.asarray(W@scores['mural']).ravel();valid=(degree>0);cov0=np.column_stack([neighbor,scores['endothelial_identity'],np.log1p(umi),scores['fibroblast'],scores['epithelial'],scores['hepatocyte']]);valid &= np.isfinite(cov0).all(1)&np.isfinite(scores['mural'])&np.isfinite(scores['hypoxia'])
  for model in (['A_prior','B_focal','C_stress'] if radius==200 else ['C_stress']):
   cv=cov0 if model=='A_prior' else np.c_[cov0,scores['mural']] if model=='B_focal' else np.c_[cov0,scores['mural'],scores['hypoxia']]
   z,sd=standard(cv[valid]);F=np.c_[np.ones(valid.sum()),z]
   if model!='A_prior' and job['cohort']=='GSE294385':F=np.c_[F,pd.get_dummies(pd.Series(regions[valid]),drop_first=True,dtype=float).to_numpy()]
   rank=np.linalg.matrix_rank(F);cond=np.linalg.cond(F);yy=Y[valid] if radius==200 else Y[valid,:1];zy,ysd=standard(yy);betas=np.linalg.lstsq(F,zy,rcond=None)[0][1] if rank==F.shape[1] and cond<1e8 else np.full(yy.shape[1],np.nan);betas[ysd<1e-12]=np.nan
   for gene,beta in zip(targets if radius==200 else ['CGNL1'],betas):rows.append(dict(cohort=job['cohort'],sample=job['sample'],patient=job['patient'],radius_um=radius,model=model,gene=gene,n_spots=int(valid.sum()),beta=float(beta),rank=int(rank),parameters=F.shape[1],condition=float(cond)))
  qc.append(dict(cohort=job['cohort'],sample=job['sample'],radius_um=radius,QC_spots=len(umi),isolated=int((degree==0).sum())))
 pd.DataFrame(rows).to_csv(O/'R1_section_effects.tsv',sep='\t',index=False);print('R1_MODEL_DONE',job['sample'],flush=True)
pd.DataFrame(qc).to_csv(O/'R1_spatial_QC.tsv',sep='\t',index=False);m=pd.DataFrame(rows);pt=m.groupby(['cohort','patient','radius_um','model','gene'],as_index=False).agg(beta=('beta','mean'),sections=('beta','count'));pt.to_csv(O/'R1_patient_effects.tsv',sep='\t',index=False)
def inference(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v);d=dict(n=n,positive=int((v>0).sum()),negative=int((v<0).sum()),mean_beta=float(v.mean()) if n else np.nan,median_beta=float(np.median(v)) if n else np.nan,P=np.nan,signflip_P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=4 and v.std()>0:
  d['P']=float(stats.ttest_1samp(v,0).pvalue);ci=stats.t.interval(.95,n-1,loc=v.mean(),scale=stats.sem(v));d.update(CI_low=float(ci[0]),CI_high=float(ci[1]));signs=np.array(list(itertools.product([-1,1],repeat=n)));d['signflip_P']=float(np.mean(np.abs(signs@v/n)>=abs(v.mean())-1e-12))
 return d
summary=[]
for key,g in pt.groupby(['cohort','radius_um','model','gene']):summary.append(dict(zip(['cohort','radius_um','model','gene'],key),**inference(g.beta)))
s=pd.DataFrame(summary);s['benchmark_FDR']=np.nan
for key,ix in s[s.radius_um.eq(200)].groupby(['cohort','model']).groups.items():
 valid=s.loc[ix].index[s.loc[ix,'P'].notna()];s.loc[valid,'benchmark_FDR']=multipletests(s.loc[valid,'P'],method='fdr_bh')[1]
s.to_csv(O/'R1_all_gene_models.tsv',sep='\t',index=False);target=s[s.gene.eq('CGNL1')].copy();target['target_family']=np.where(target.radius_um.eq(200),'R1_target_models3','R1_target_radius_sensitivity2');target['target_FDR']=np.nan
for key,ix in target.groupby(['cohort','target_family']).groups.items():
 valid=target.loc[ix].index[target.loc[ix,'P'].notna()];target.loc[valid,'target_FDR']=multipletests(target.loc[valid,'P'],method='fdr_bh')[1]
target.to_csv(O/'R1_target_summary.tsv',sep='\t',index=False);specific=[];ranks=[]
for (cohort,model),g in pt[pt.radius_um.eq(200)].groupby(['cohort','model']):
 pivot=g.pivot(index='patient',columns='gene',values='beta');delta=pivot.CGNL1-pivot.drop(columns='CGNL1').median(axis=1);specific.append(dict(cohort=cohort,model=model,contrast='target_minus_median_matched_control',**inference(delta)))
 means=pivot.mean();v=means.drop('CGNL1');ranks.append(dict(cohort=cohort,model=model,target_mean=float(means.CGNL1),controls_median=float(v.median()),signed_percentile=float(np.mean(v<means.CGNL1)),absolute_percentile=float(np.mean(np.abs(v)<abs(means.CGNL1))),n_controls=len(v),rank_is_not_biological_P=True))
t=pd.DataFrame(specific);t['FDR']=np.nan
for cohort,ix in t.groupby('cohort').groups.items():
 valid=t.loc[ix].index[t.loc[ix,'P'].notna()];t.loc[valid,'FDR']=multipletests(t.loc[valid,'P'],method='fdr_bh')[1]
t.to_csv(O/'R1_specificity_contrast.tsv',sep='\t',index=False);pd.DataFrame(ranks).to_csv(O/'R1_control_rank.tsv',sep='\t',index=False)
(O/'SPATIAL_SPECIFICITY_COMPLETE.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(jobs),'controls':len(chosen),'cohorts':['GSE294385','GSE217414'],'no_spot_inferential_tests':True,'original_results_modified':False},indent=2));print('SPATIAL_SPECIFICITY_COMPLETE',flush=True)
