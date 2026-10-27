from pathlib import Path
import json,gzip,datetime,hashlib,itertools
import numpy as np,pandas as pd
from scipy import sparse,io,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];B=R.parent;D=R/'results';C=json.loads((R/'SPECIFICATION.json').read_text());seed=20261002;rng=np.random.default_rng(seed)
identity=B/'identity_review_target_blind_20261002_v1';accepted=pd.read_csv(R/'accepted_review/accepted_cluster_annotation.tsv',sep='\t',dtype=str).set_index(['original_dataset','original_cluster']);control=pd.read_csv(B/'formal_v1_20261002/05_detection_controls/discovery_selected_controls.tsv',sep='\t').iloc[:,0].tolist();assert len(control)==20
ann={};cache=[];coverage=[];profiles=[];matchrows=[];controlrows=[]
for ds in C['endothelial_clusters']:
 t=pd.read_csv(identity/'results'/ds/'cell_identity_evidence.tsv.gz',sep='\t',dtype={'cluster05':str});t['accepted_major_lineage']=[accepted.loc[(ds,c),'accepted_major_lineage'] for c in t.cluster05];t['accepted_state']=[accepted.loc[(ds,c),'accepted_state'] for c in t.cluster05];t['mito_QC_legacy']=t.mito_QC;t['mito_QC_current']=np.nan if ds=='GSE245552' else t.mito_QC;t['mito_QC_status']='unavailable_feature_matching' if ds=='GSE245552' else 'measurable';assert t.accepted_major_lineage.notna().all();t.to_csv(D/(ds+'_accepted_annotation.tsv.gz'),sep='\t',index=False);ann[ds]=t
I=B/'spatial_v3/background_validation/run_20260924T052010Z';jobs=[]
for s in json.loads((I/'SAMPLES.json').read_text()):
 if s['dataset']=='GSE285991':jobs.append((s['dataset'],s['sample'],s['patient_label_provisional'],[I/x for x in s['files']]))
for s in json.loads((B/'cohort_verification_20260924/GSE245552_LM_SAMPLES.json').read_text()):jobs.append(('GSE245552',s['sample'],'GSE245552:'+s['donor'],[B/'cohort_verification_20260924/GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]))
for ds,sample,donor,paths in jobs:
 fs={k:next(p for p in paths if k in p.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(fs['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(fs['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
 with gzip.open(fs['matrix'],'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 X.eliminate_zeros();assert X.shape==(len(sy),len(bc)) and np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data));t=ann[ds];g=t[(t['sample']==sample)&(t.cluster05==C['endothelial_clusters'][ds])].copy();ix=pd.Index(bc).get_indexer(g.barcode);assert (ix>=0).all();Y=X[:,ix].tocsr();lib=np.asarray(Y.sum(axis=0)).ravel();ng=Y.getnnz(axis=0);assert np.allclose(lib,g.UMI_QC.to_numpy());mito_idx=np.char.startswith(sy,'MT-');mt=np.asarray(Y[mito_idx].sum(axis=0)).ravel()/np.maximum(lib,1)*100 if mito_idx.any() else np.full(len(ix),np.nan)
 name=ds+'_'+sample;sparse.save_npz(R/'cache'/(name+'.npz'),Y);pd.DataFrame({'gene':sy}).to_csv(R/'cache'/(name+'_genes.tsv'),sep='\t',index=False);g['mito_retest']=mt;g.to_csv(R/'cache'/(name+'_cells.tsv.gz'),sep='\t',index=False)
 cache.append(dict(dataset=ds,sample=sample,donor=donor,path=str(R/'cache'/(name+'.npz')),genes=str(R/'cache'/(name+'_genes.tsv')),cells=str(R/'cache'/(name+'_cells.tsv.gz')),selected_cells=len(g),mt_feature_count=int(mito_idx.sum()),source_matrix=str(fs['matrix'])))
 print('EXTRACTED',ds,sample,len(g),flush=True)
(R/'CACHE_INDEX.json').write_text(json.dumps(cache,indent=2))
def matched(y,genes,lib,ng,exposure='CGNL1',caliper=.25,minpairs=3):
 idx={g:i for i,g in enumerate(genes)}
 if exposure not in idx:return dict(pairs=0,effect=np.nan,eligible=False,reason='gene_unavailable')
 v=y[idx[exposure]].toarray().ravel();pos=np.flatnonzero(v>0);neg=np.flatnonzero(v==0);score=np.log1p(y[[idx[g] for g in C['panel']]].toarray()/lib*1e4).mean(axis=0)
 if caliper is None:
  ok=len(pos)>=minpairs and len(neg)>=minpairs;return dict(pairs=min(len(pos),len(neg)),effect=float(score[pos].mean()-score[neg].mean()) if ok else np.nan,eligible=ok,reason='' if ok else 'insufficient_groups')
 feat=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for i in sorted(pos,key=lambda j:(v[j],j)):
  cand=np.array([j for j in neg if j not in used],int)
  if not len(cand):break
  delta=np.abs(feat[cand]-feat[i]);ok=(delta<=caliper).all(axis=1)
  if ok.any():
   cand=cand[ok];delta=delta[ok];j=cand[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((i,j))
 if len(pairs)<minpairs:return dict(pairs=len(pairs),effect=np.nan,eligible=False,reason='insufficient_matched_pairs')
 a,b=np.array(pairs).T;return dict(pairs=len(pairs),effect=float((score[a]-score[b]).mean()),eligible=True,reason='',mean_abs_logUMI_difference=float(np.abs(feat[a,0]-feat[b,0]).mean()),mean_abs_logGene_difference=float(np.abs(feat[a,1]-feat[b,1]).mean()))
for rec in cache:
 Y=sparse.load_npz(rec['path']);genes=pd.read_csv(rec['genes'],sep='\t').gene.to_numpy(dtype=str);g=pd.read_csv(rec['cells'],sep='\t');lib=np.asarray(Y.sum(axis=0)).ravel();ng=Y.getnnz(axis=0);mt=g.mito_retest.to_numpy();idx={gene:i for i,gene in enumerate(genes)}
 for gene in ['CGNL1']+C['panel']:assert np.sum(genes==gene)==1
 for s in C['continuous_stages']:
  k=(lib>=s['umi'])&(ng>=s['genes'])&(np.isnan(mt)|(mt<=s['mt']))
  if s['branch']=='exclude_reference_immune':k &= ~g.reference_agreement.eq('both_immune').to_numpy()
  elif s['branch']=='both_references_endo':k &= g.reference_agreement.eq('both_endothelial').to_numpy()
  n=int(k.sum());base=dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],stage=s['stage'],cells=n,eligible=n>=s['mincells'],mito_QC_available=rec['mt_feature_count']>0);coverage.append(base)
  if n>=s['mincells']:
   small=Y[[idx[x] for x in ['CGNL1']+C['panel']]][:,k].toarray();v=np.log1p(small.sum(axis=1)/lib[k].sum()*1e4);profiles.append(base|dict(CGNL1=v[0],junction=float(v[1:].mean()),mean_logUMI=float(np.log1p(lib[k]).mean()),mean_logGenes=float(np.log1p(ng[k]).mean()),CGNL1_counts=float(small[0].sum()),CGNL1_positive_cells=int((small[0]>0).sum()),**{gene+'_pb':float(v[i]) for i,gene in enumerate(['CGNL1']+C['panel'])}))
 for s in C['match_stages']:
  res=matched(Y,genes,lib,ng,caliper=s['caliper'],minpairs=s['minpairs']) if len(g)>=s['mincells'] else dict(pairs=0,effect=np.nan,eligible=False,reason='insufficient_cells')
  matchrows.append(dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],stage=s['stage'],cells=len(g),**res))
 if len(g)>=10:
  for gene in ['CGNL1']+control:controlrows.append(dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],gene=gene,**matched(Y,genes,lib,ng,exposure=gene)))
coverage=pd.DataFrame(coverage);coverage.to_csv(D/'lesion_eligibility.tsv',sep='\t',index=False);pr=pd.DataFrame(profiles);pr.to_csv(D/'lesion_continuous_profiles.tsv',sep='\t',index=False);pt=pr.groupby(['stage','dataset','donor'],as_index=False).mean(numeric_only=True);pt.to_csv(D/'patient_continuous_profiles.tsv',sep='\t',index=False)
mr=pd.DataFrame(matchrows);mr.to_csv(D/'matched_lesion_results.tsv',sep='\t',index=False);mp=mr[mr.eligible].groupby(['stage','dataset','donor'],as_index=False).effect.mean();mp.to_csv(D/'matched_patient_effects.tsv',sep='\t',index=False)
def select(v,scope):return v if scope=='pooled' else v[v.dataset==scope]
def fit(v,model,cov='nonrobust'):
 result=dict(n=len(v),n_studies=v.dataset.nunique(),beta=np.nan,slope=np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,status='not_estimable')
 if len(v)<4:return result|dict(reason='n_below4')
 x=v.CGNL1.to_numpy();y=v.junction.to_numpy()
 if np.std(x)<1e-12 or np.std(y)<1e-12:return result|dict(reason='zero_variation')
 F=pd.DataFrame({'CGNL1':x},index=v.index)
 if model!='unadjusted':F=pd.concat([F,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
 if model=='study_depth_adjusted':F['depth']=v.mean_logUMI
 F=sm.add_constant(F,has_constant='add')
 if np.linalg.matrix_rank(F.to_numpy())<F.shape[1] or len(v)-F.shape[1]<3:return result|dict(reason='rank_or_residual_df')
 m=sm.OLS(y,F).fit(cov_type=cov,use_t=True);ci=m.conf_int().loc['CGNL1'];return result|dict(beta=float(m.params['CGNL1'])*x.std(ddof=1)/y.std(ddof=1),slope=float(m.params['CGNL1']),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(m.pvalues['CGNL1']),status='estimated',reason='',df=float(m.df_resid))
models=[];scopes=['pooled','GSE245552','GSE285991']
for s in C['continuous_stages']:
 for scope in scopes:
  v=select(pt[pt.stage==s['stage']],scope)
  for model in C['continuous_models']:
   primary=s['stage']=='cluster_primary_n10' and model==('study_adjusted' if scope=='pooled' else 'unadjusted');models.append(dict(stage=s['stage'],scope=scope,model=model,is_primary=primary,family='A05_primary_3' if primary else 'A05_sensitivity_'+scope,**fit(v,model)))
def fdrtable(data,name):
 t=pd.DataFrame(data);t['FDR']=np.nan
 for fam,ids in t.groupby('family').groups.items():
  valid=t.loc[ids].index[t.loc[ids,'P'].notna()]
  if len(valid):t.loc[valid,'FDR']=multipletests(t.loc[valid,'P'],method='fdr_bh')[1]
 t['nominal_P005']=t.P<.05;t.to_csv(D/name,sep='\t',index=False);return t
res=fdrtable(models,'continuous_models.tsv');rob=[];interaction=[]
for scope in scopes:
 v=select(pt[pt.stage=='cluster_primary_n10'],scope).reset_index(drop=True);model='study_adjusted' if scope=='pooled' else 'unadjusted';rr=fit(v,model,cov='HC3');rob.append(dict(scope=scope,method='HC3',family='A05_robustness',**rr))
 if rr['status']=='estimated':
  Z=sm.add_constant(pd.get_dummies(v.dataset,drop_first=True,dtype=float),has_constant='add').to_numpy();x=v.CGNL1.to_numpy();y=v.junction.to_numpy();xr=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];yr=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0];obs=abs(xr@yr);extreme=0;groups=[np.flatnonzero(v.dataset.to_numpy()==ds) for ds in v.dataset.unique()]
  for i in range(4999):
   yp=yr.copy()
   for ids in groups:yp[ids]=rng.permutation(yr[ids])
   extreme += abs(xr@yp)>=obs-1e-12
  rr=fit(v,model);rr['P']=(extreme+1)/5000;rob.append(dict(scope=scope,method='within_study_permutation_CI_OLS',family='A05_robustness',**rr))
fdrtable(rob,'continuous_robustness.tsv')
v=pt[pt.stage=='cluster_primary_n10'].copy();dum=pd.get_dummies(v.dataset,drop_first=True,dtype=float);F=sm.add_constant(pd.concat([v[['CGNL1']],dum],axis=1),has_constant='add');G=F.copy()
for col in dum:G['interaction_'+col]=v.CGNL1*dum[col]
if np.linalg.matrix_rank(G.to_numpy())==G.shape[1] and len(v)-G.shape[1]>=3:
 m0=sm.OLS(v.junction,F).fit();m1=sm.OLS(v.junction,G).fit();f,p,df=m1.compare_f_test(m0);pd.DataFrame([dict(n=len(v),F=f,P=p,df_difference=df,residual_df=m1.df_resid,family='single_prespecified_study_interaction',FDR=np.nan)]).to_csv(D/'study_slope_interaction.tsv',sep='\t',index=False)
matched_models=[]
for s in C['match_stages']:
 for scope in scopes:
  v=select(mp[mp.stage==s['stage']],scope);x=v.effect.to_numpy();n=len(x);primary=s['stage']=='matched_primary_n10_pairs3';row=dict(stage=s['stage'],scope=scope,n=n,positive=int((x>0).sum()),effect=float(x.mean()) if n else np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,signflip_P=np.nan,is_primary=primary,family='A06_primary' if primary else 'A06_sensitivity_'+scope)
  if n>=4 and np.std(x,ddof=1)>1e-12:
   ci=stats.t.interval(.95,n-1,loc=x.mean(),scale=stats.sem(x));row.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=float(ci[0]),CI_high=float(ci[1]))
   if n<=18:
    observed=abs(x.mean());extreme=sum(abs(np.dot(signs,x)/n)>=observed-1e-12 for signs in itertools.product([-1,1],repeat=n));row['signflip_P']=extreme/(2**n)
  matched_models.append(row)
ms=fdrtable(matched_models,'matched_models.tsv')
ct=pd.DataFrame(controlrows);ct.to_csv(D/'technical_control_lesion_results.tsv',sep='\t',index=False);cp=ct[ct.eligible].groupby(['dataset','donor','gene'],as_index=False).effect.mean();cp.to_csv(D/'technical_control_patient_results.tsv',sep='\t',index=False);target=cp[cp.gene=='CGNL1'].set_index('donor').effect;compare=[]
for gene in control:
 vv=cp[cp.gene==gene].set_index('donor').effect;common=target.index.intersection(vv.index)
 if len(common):compare.append(dict(gene=gene,n_common=len(common),CGNL1_mean=float(target[common].mean()),control_mean=float(vv[common].mean()),difference=float((target[common]-vv[common]).mean()),positive_patient_differences=int((target[common]>vv[common]).sum()),scope='descriptive_not_a_null_or_causal_test'))
pd.DataFrame(compare).to_csv(D/'technical_control_common_patient_comparison.tsv',sep='\t',index=False)
# Inspect input identities before any conclusion about significance.
primary=pt[pt.stage=='cluster_primary_n10'];assert not primary[['dataset','donor']].duplicated().any();summary={'exit_code':0,'accepted_clusters':C['endothelial_clusters'],'primary_patient_n':len(primary),'patient_counts_by_dataset':primary.dataset.value_counts().to_dict(),'continuous_primary':res[res.is_primary].to_dict('records'),'matched_primary':ms[ms.is_primary].to_dict('records'),'mito_QC_GSE245552':'unavailable, not zero; no MT-dependent exclusion applied','finished':datetime.datetime.now().astimezone().isoformat()};(R/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
