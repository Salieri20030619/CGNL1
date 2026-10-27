from pathlib import Path
import json,gzip,datetime,hashlib,re
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse,io,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];B=R.parent;C=json.loads((R/'SPECIFICATION.json').read_text());D=R/'results';genes5=['CGNL1']+C['panel'];cache=[];checks=[];profiles=[];matchrows=[];coverage=[]
labels=list(C['identity_markers']);ej=labels.index('endothelial')
def emit(ds,donor,sample,X,sy,bc,infer=True,author=None,qfields=None,adjusted=False):
 X=X.tocsr().astype(np.float64);X.eliminate_zeros();sy=np.asarray(sy,dtype=str);assert X.shape==(len(sy),len(bc));assert len(set(bc))==len(bc);assert np.isfinite(X.data).all() and (X.data>=0).all()
 if not adjusted:assert np.allclose(X.data,np.round(X.data))
 lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100
 scores=[];hit0=[];hit01=[];hit1=[]
 for gs in C['identity_markers'].values():
  zz=X[np.isin(sy,gs)].toarray();assert len(zz)>0;scores.append(np.log1p(zz/np.maximum(lib,1)*1e4).mean(axis=0));hit0.append((zz>0).sum(axis=0));hit01.append((zz>=.1).sum(axis=0));hit1.append((zz>=1).sum(axis=0))
 scores=np.array(scores).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2]
 indices=[];measurement_ok=True
 for g in genes5:
  ix=np.flatnonzero(sy==g);measurement_ok &= len(ix)==1;indices.append(int(ix[0]) if len(ix)==1 else -1)
 qlib,qng,qmt=(lib,ng,mt) if qfields is None else qfields
 if adjusted:base=np.asarray(author,bool)&(lib>0)
 else:base=(best==ej)&(np.array(hit0)[ej]>=1)&(lib>0)
 checks.append(dict(dataset=ds,donor=donor,sample=sample,total_cells=len(bc),candidate_cells=int(base.sum()),infer_patient=infer,scale='decontX_adjusted' if adjusted else 'raw_counts'))
 z=dict(dataset=ds,donor=donor,sample=sample,infer=infer,adjusted=adjusted,counts=np.array([X[i,base].toarray().ravel() if i>=0 else np.full(int(base.sum()),np.nan) for i in indices]),lib=lib[base],ng=ng[base],qlib=np.asarray(qlib)[base],qng=np.asarray(qng)[base],qmt=np.asarray(qmt)[base],hit=np.array(hit0)[ej][base],hit01=np.array(hit01)[ej][base],hit1=np.array(hit1)[ej][base],margin=margin[base],best_endo=(best[base]==ej),barcode=np.asarray(bc,dtype=str)[base])
 f=R/'cache'/(ds+'_'+sample+'.npz');np.savez_compressed(f,**{k:v for k,v in z.items() if k not in ['dataset','donor','sample','infer','adjusted']});cache.append({k:z[k] for k in ['dataset','donor','sample','infer','adjusted']}|{'path':str(f)})
 print('CACHE',ds,sample,len(z['lib']),flush=True)
def load(paths):
 fs={k:next(p for p in paths if k in p.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(fs['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(fs['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
 with gzip.open(fs['matrix'],'rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsr()
 return X,sy,bc
# Full original cohorts, including previously excluded libraries.
X=sparse.load_npz(B/'deep_v2/results/cnv/counts_gene_by_cell.npz');sy=pd.read_csv(B/'deep_v2/results/cnv/gene_names.tsv',sep='\t').gene.to_numpy(dtype=str);m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);assert X.shape==(len(sy),len(m))
for pt,g in m[m.organs=='LCT'].groupby('patients'):
 ix=m.index.get_indexer(g.index);emit('GSE225857','GSE225857:'+pt,pt,X[:,ix],sy,m.index.to_numpy()[ix])
del X
I=B/'spatial_v3/background_validation/run_20260924T052010Z'
for s in json.loads((I/'SAMPLES.json').read_text()):
 if s['dataset']=='GSE231559' and not (s['title'].startswith('L') and not s['title'].endswith('N')):continue
 X,sy,bc=load([I/p for p in s['files']]);emit(s['dataset'],s['patient_label_provisional'],s['sample'],X,sy,bc)
old=B/'cohort_verification_20260924'
for s in json.loads((old/'GSE245552_LM_SAMPLES.json').read_text()):
 X,sy,bc=load([old/'GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]);emit('GSE245552','GSE245552:'+s['donor'],s['sample'],X,sy,bc)
p=old/'GSE178318';sy=pd.read_csv(p/'GSE178318_genes.tsv.gz',sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(p/'GSE178318_barcodes.tsv.gz',sep='\t',header=None)[0].to_numpy(dtype=str)
with gzip.open(p/'GSE178318_matrix.mtx.gz','rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsr()
tags=np.array([x.split('_',1)[1] for x in bc])
for tag in sorted(set(tags)):
 if tag.endswith('_LM'):
  ix=np.flatnonzero(tags==tag);emit('GSE178318','GSE178318:'+tag.split('_')[0],tag,X[:,ix],sy,bc[ix])
del X
for f in sorted((B/'data_prepare_20260929/wave01_core/GSE234804').glob('*_LM*')):
 with h5py.File(f) as h:
  sy=h['assays/RNA/features'].asstr()[:];bc=h['cell.names'].asstr()[:];arr=h['assays/RNA/counts'];X=sparse.vstack([sparse.csr_matrix(arr[i:i+256]) for i in range(0,len(bc),256)],format='csr')
 emit('GSE234804','UNRESOLVED:'+f.name.split('_')[0],f.name.split('_')[0],X.T.tocsr(),sy,bc,infer=False)
a=ad.read_h5ad(B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad')
for pt in a.obs.Patient.astype(str).unique():
 v=a[a.obs.Patient.astype(str)==pt];q=v.obs;emit('multiome','multiome:'+pt,pt,v.X.T.tocsr(),v.var_names.to_numpy(dtype=str),v.obs_names.to_numpy(dtype=str),author=q.Cell_type.astype(str).eq('Endothelial').to_numpy(),qfields=(q.nCount_RNA.to_numpy(),q.nFeature_RNA.to_numpy(),q['percent.mt'].to_numpy()),adjusted=True)
(R/'CACHE_INDEX.json').write_text(json.dumps(cache,indent=2));pd.DataFrame(checks).to_csv(D/'all_library_inventory.tsv',sep='\t',index=False)
# Freeze dimensions in cache, then run the complete predeclared ladder.
def mask(z,s,adj):
 k=(z['qlib']>=s['umi'])&(z['qng']>=s['genes'])&(z['qmt']<=s['mt'])
 if adj:
  if s['independent_best']:k&=z['best_endo']
  hits=z['hit01'] if s['hit_threshold']==.1 else z['hit1'];k&=hits>=s['hits']
 else:k&=z['hit']>=s['hits']
 return k&(z['margin']>=s['margin'] if not adj or s['independent_best'] else True)
def match_effect(z,k,caliper,minpairs):
 c=z['counts'][:,k];lib=z['lib'][k];ng=z['ng'][k];pos=np.flatnonzero(c[0]>0);neg=np.flatnonzero(c[0]==0);v=np.log1p(c[1:]/lib*1e4).mean(axis=0)
 if caliper is None:
  eligible=len(pos)>=minpairs and len(neg)>=minpairs;return (float(v[pos].mean()-v[neg].mean()) if eligible else np.nan,min(len(pos),len(neg)),eligible)
 feat=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for i in sorted(pos,key=lambda j:(c[0,j],j)):
  cand=np.array([j for j in neg if j not in used],int)
  if not len(cand):break
  dist=np.abs(feat[cand]-feat[i]);ok=(dist<=caliper).all(axis=1)
  if ok.any():
   cand=cand[ok];dist=dist[ok];j=cand[np.argmin((dist**2).sum(axis=1))];used.add(j);pairs.append((i,j))
 if len(pairs)<minpairs:return np.nan,len(pairs),False
 a,b=np.array(pairs).T;return float((v[a]-v[b]).mean()),len(pairs),True
for v in cache:
 with np.load(v['path']) as h:z={k:h[k] for k in h.files}
 stages=C['multiome_stages'] if v['adjusted'] else C['raw_stages'];saved={}
 for s in stages:
  k=mask(z,s,v['adjusted']);n=int(k.sum());saved[s['stage']]=k;base=dict(dataset=v['dataset'],donor=v['donor'],sample=v['sample'],infer=v['infer'],branch='multiome' if v['adjusted'] else 'raw',stage=s['stage'],cells=n,eligible=n>=s['mincells'] and np.isfinite(z['counts']).all(),measurement_ok=bool(np.isfinite(z['counts']).all()))
  coverage.append(base)
  if base['eligible']:
   c=z['counts'][:,k];lib=z['lib'][k];pb=np.log1p(c.sum(axis=1)/lib.sum()*1e4);profiles.append(base|dict(CGNL1=pb[0],junction=float(pb[1:].mean()),mean_logUMI=float(np.log1p(z['qlib'][k]).mean())))
 if not v['adjusted']:
  for s in C['match_stages']:
   cs=next(x for x in C['raw_stages'] if x['stage']==s['cell_stage']);k=saved[s['cell_stage']];n=int(k.sum())
   effect,pairs,eligible=match_effect(z,k,s['caliper'],s['minpairs']) if n>=cs['mincells'] and np.isfinite(z['counts']).all() else (np.nan,0,False)
   matchrows.append(dict(dataset=v['dataset'],donor=v['donor'],sample=v['sample'],infer=v['infer'],stage=s['stage'],cell_stage=s['cell_stage'],cells=n,pairs=pairs,effect=effect,eligible=eligible))
pd.DataFrame(coverage).to_csv(D/'library_eligibility_all_stages.tsv',sep='\t',index=False);pr=pd.DataFrame(profiles);pr.to_csv(D/'lesion_profiles_all_stages.tsv',sep='\t',index=False)
pat=pr.groupby(['branch','stage','dataset','donor','infer'],as_index=False)[['CGNL1','junction','mean_logUMI','cells']].mean();pat.to_csv(D/'patient_profiles_all_stages.tsv',sep='\t',index=False)
mt=pd.DataFrame(matchrows);mt.to_csv(D/'matching_library_all_stages.tsv',sep='\t',index=False);mp=mt[mt.eligible&mt.infer].groupby(['stage','dataset','donor'],as_index=False).effect.mean();mp.to_csv(D/'matching_patient_all_stages.tsv',sep='\t',index=False)
def fit(v,model):
 result=dict(n_patients=len(v),n_studies=v.dataset.nunique(),beta=np.nan,slope=np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,status='not_estimable')
 if len(v)<4:return result|{'reason':'fewer_than4_patients'}
 x=v.CGNL1.to_numpy();y=v.junction.to_numpy()
 if np.std(x)<1e-12 or np.std(y)<1e-12:return result|{'reason':'zero_variation'}
 design=pd.DataFrame({'CGNL1':x},index=v.index)
 if model!='M0_unadjusted':design=pd.concat([design,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
 if model=='M2_study_depth':design['depth']=v.mean_logUMI
 design=sm.add_constant(design,has_constant='add')
 if np.linalg.matrix_rank(design.to_numpy())<design.shape[1]:return result|{'reason':'rank_deficient'}
 if len(v)-design.shape[1]<3:return result|{'reason':'residual_df_below3'}
 mod=sm.OLS(y,design).fit();ci=mod.conf_int().loc['CGNL1'];slope=float(mod.params['CGNL1']);return result|dict(beta=slope*np.std(x,ddof=1)/np.std(y,ddof=1),slope=slope,CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(mod.pvalues['CGNL1']),status='estimated',reason='',df=float(mod.df_resid))
def subset(v,scope):
 if scope=='all':return v
 if scope=='external':return v[v.dataset!='GSE225857']
 return v[v.dataset==scope]
models=[];decomp=[];transitions=[]
for branch,stages,scopes in [('raw',C['raw_stages'],['all','external','GSE225857','GSE231559','GSE285991','GSE245552','GSE178318']),('multiome',C['multiome_stages'],['all'])]:
 a=pat[(pat.branch==branch)&pat.infer]
 for scope in scopes:
  prev=None
  for s in stages:
   v=subset(a[a.stage==s['stage']],scope)
   for model in C['models']:
    models.append(dict(branch=branch,scope=scope,stage=s['stage'],model=model,family=branch+'_continuous_'+scope,**fit(v,model)))
    if prev is not None:
     u=subset(a[a.stage==prev],scope);common=set(u.donor)&set(v.donor);uu=u[u.donor.isin(common)];vv=v[v.donor.isin(common)];f0=fit(uu,model);f1=fit(vv,model)
     decomp.append(dict(branch=branch,scope=scope,previous=prev,current=s['stage'],model=model,n_common=len(common),lost_patients=';'.join(sorted(set(u.donor)-set(v.donor))),added_patients=';'.join(sorted(set(v.donor)-set(u.donor))),previous_common_slope=f0['slope'],current_common_slope=f1['slope'],previous_common_beta=f0['beta'],current_common_beta=f1['beta'],note='common patient comparison; remaining change includes cell selection and qualifying-lesion composition'))
   prev=s['stage']
res=pd.DataFrame(models);res['FDR']=np.nan
for family,ix in res.groupby('family').groups.items():
 valid=res.loc[ix,'P'].notna();ids=res.loc[ix].index[valid]
 if len(ids):res.loc[ids,'FDR']=multipletests(res.loc[ids,'P'],method='fdr_bh')[1]
res['nominal']=res.P<.05;res.to_csv(D/'continuous_models_all_stages.tsv',sep='\t',index=False);pd.DataFrame(decomp).to_csv(D/'common_patient_decomposition.tsv',sep='\t',index=False)
for (branch,scope,model),v in res.groupby(['branch','scope','model'],sort=False):
 prev=None
 for x in v.to_dict('records'):
  if prev is not None and np.isfinite(prev['P']) and np.isfinite(x['P']) and prev['nominal']!=x['nominal']:transitions.append(dict(branch=branch,scope=scope,model=model,from_stage=prev['stage'],to_stage=x['stage'],transition='gain' if x['nominal'] else 'loss',P_before=prev['P'],P_after=x['P'],n_before=prev['n_patients'],n_after=x['n_patients']))
  prev=x
pd.DataFrame(transitions,columns=['branch','scope','model','from_stage','to_stage','transition','P_before','P_after','n_before','n_after']).to_csv(D/'continuous_nominal_transitions.tsv',sep='\t',index=False)
ms=[]
for scope in ['all','external','GSE225857']:
 for s in C['match_stages']:
  v=subset(mp[mp.stage==s['stage']],scope);x=v.effect.to_numpy();n=len(x);row=dict(scope=scope,stage=s['stage'],n_patients=n,positive=int((x>0).sum()),effect=float(x.mean()) if n else np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,family='within_patient_'+scope)
  if n>=4 and np.std(x,ddof=1)>1e-12:
   se=stats.sem(x);ci=stats.t.interval(.95,n-1,loc=x.mean(),scale=se);row.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=float(ci[0]),CI_high=float(ci[1]))
  ms.append(row)
ms=pd.DataFrame(ms);ms['FDR']=np.nan
for fam,ix in ms.groupby('family').groups.items():
 ids=ms.loc[ix].index[ms.loc[ix,'P'].notna()]
 if len(ids):ms.loc[ids,'FDR']=multipletests(ms.loc[ids,'P'],method='fdr_bh')[1]
ms['nominal']=ms.P<.05;ms.to_csv(D/'within_patient_models_all_stages.tsv',sep='\t',index=False)
# Endpoint reconciliation with the existing strict branch, not a selective rerun.
strict=pat[(pat.branch=='raw')&(pat.stage=='N2_cells50')&pat.infer];prior=pd.read_csv(B/'formal_v1_20261002/04_patient_state/patient_pseudobulk_profiles.tsv',sep='\t');j=prior.merge(strict,on=['dataset','donor'],suffixes=('_old','_new'));j['CGNL1_delta']=j.CGNL1-j.CGNL1_pb;j['junction_delta']=j.junction-j.junction_pb;j.to_csv(D/'strict_endpoint_reconciliation.tsv',sep='\t',index=False)
summary={'exit_code':0,'full_stage_grid_completed':True,'raw_libraries':sum(not x['adjusted'] for x in cache),'multiome_patients':sum(x['adjusted'] for x in cache),'continuous_model_rows':len(res),'within_patient_model_rows':len(ms),'strict_final_n':len(strict),'strict_prior_n':len(prior),'strict_overlap_n':len(j),'strict_max_target_difference':float(j.CGNL1_delta.abs().max()),'strict_max_score_difference':float(j.junction_delta.abs().max()),'continuous_threshold_crossings':len(transitions),'finished':datetime.datetime.now().astimezone().isoformat()}
(R/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
