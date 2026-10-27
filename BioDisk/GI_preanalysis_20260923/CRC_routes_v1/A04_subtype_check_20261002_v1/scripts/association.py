
from pathlib import Path
import json,hashlib,datetime,itertools
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'A04_subtype_check_20261002_v1';D=R/'results';D.mkdir(exist_ok=True)
freeze=json.loads((R/'IDENTITY_FREEZE.json').read_text());assert all(hashlib.sha256((R/'identity'/f).read_bytes()).hexdigest()==h for f,h in freeze['hashes'].items())
cache=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text());ann={ds:pd.read_csv(R/'identity'/(ds+'_cell_subtype.tsv.gz'),sep='\t').set_index('cell_id') for ds in ['GSE245552','GSE285991']}
panel=['CLDN5','ESAM','TJP1','OCLN'];types=['sinusoidal_like','other_blood_endo_like','lymphatic_like'];profiles=[];subprofiles=[];matches=[];detection=[]
def pb(Y,lib,ng,k):
 v=np.log1p(Y[:,k].sum(axis=1)/lib[k].sum()*1e4)
 return dict(CGNL1=float(v[0]),junction=float(v[1:].mean()),mean_logUMI=float(np.log1p(lib[k]).mean()),cells=int(k.sum()),counts=int(Y[0,k].sum()),positive_cells=int((Y[0,k]>0).sum()),**{g:float(v[i+1]) for i,g in enumerate(panel)})
def match(Y,lib,ng,labels,restrict,minpairs):
 v=Y[0];score=np.log1p(Y[1:]/lib*1e4).mean(axis=0);pos=np.flatnonzero(v>0);neg=np.flatnonzero(v==0);feat=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for i in sorted(pos,key=lambda j:(v[j],j)):
  cand=np.array([j for j in neg if j not in used and (not restrict or labels[j]==labels[i])],int)
  if not len(cand):continue
  delta=np.abs(feat[cand]-feat[i]);ok=(delta<=.25).all(axis=1)
  if ok.any():
   cand=cand[ok];delta=delta[ok];j=cand[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((i,j))
 ok=len(pairs)>=minpairs
 if not ok:return dict(pairs=len(pairs),eligible=False,effect=np.nan)
 a,b=np.array(pairs).T;return dict(pairs=len(pairs),eligible=True,effect=float((score[a]-score[b]).mean()),mean_abs_logUMI_difference=float(np.abs(feat[a,0]-feat[b,0]).mean()),mean_abs_logGenes_difference=float(np.abs(feat[a,1]-feat[b,1]).mean()),subtype_mismatch_pairs=int(np.sum(labels[a]!=labels[b])))
for rec in cache:
 if rec['selected_cells']==0:continue
 g=pd.read_csv(rec['cells'],sep='\t');o=ann[rec['dataset']].loc[g.cell_id];assert (o.barcode.to_numpy()==g.barcode.to_numpy()).all();full=sparse.load_npz(rec['path']);lib=np.asarray(full.sum(axis=0)).ravel();ng=full.getnnz(axis=0);genes=pd.read_csv(rec['genes'],sep='\t').gene;idx=[int(np.flatnonzero(genes.eq(x))[0]) for x in ['CGNL1']+panel];Y=full[idx].toarray();labels=o.coarse_subtype.to_numpy();base=dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'])
 for branch in ['all','prior_reference_supported']:
  k=np.ones(len(g),bool) if branch=='all' else g.reference_agreement.eq('both_endothelial').to_numpy();nk=int(k.sum());bb=base|dict(branch=branch)
  if nk>=10:
   profiles.append(bb|pb(Y,lib,ng,k)|dict(sinusoidal_fraction=float(np.mean(labels[k]=='sinusoidal_like')),lymphatic_fraction=float(np.mean(labels[k]=='lymphatic_like'))))
   detection.append(bb|dict(cells=nk,target_counts=int(Y[0,k].sum()),target_positive=int((Y[0,k]>0).sum()),one_molecule_positive=int((Y[0,k]==1).sum()),max_target_count=int(Y[0,k].max())))
   for constraint in ['depth_only','same_coarse_subtype']:
    for minpairs in [3,10]:matches.append(bb|dict(constraint=constraint,minpairs=minpairs,cells=nk,**match(Y[:,k],lib[k],ng[k],labels[k],constraint!='depth_only',minpairs)))
  for mincells in ([3,5,10] if branch=='all' else [5]):
   for label in types:
    kk=k&(labels==label)
    if kk.sum()>=mincells:subprofiles.append(bb|dict(mincells=mincells,subtype=label,**pb(Y,lib,ng,kk)))
def patients(t,keys):
 metrics=['CGNL1','junction','mean_logUMI']+panel+[c for c in ['sinusoidal_fraction','lymphatic_fraction'] if c in t]
 ag={c:'mean' for c in metrics}|{'cells':'sum','counts':'sum','positive_cells':'sum','sample':'count'}
 return t.groupby(keys,as_index=False).agg(ag).rename(columns={'sample':'qualifying_lesions'})
pr=pd.DataFrame(profiles);pr.to_csv(D/'lesion_profiles_composition.tsv',sep='\t',index=False);pt=patients(pr,['branch','dataset','donor']);pt.to_csv(D/'patient_profiles_composition.tsv',sep='\t',index=False)
sp=pd.DataFrame(subprofiles);sp.to_csv(D/'lesion_subtype_profiles.tsv',sep='\t',index=False);st=patients(sp,['branch','mincells','subtype','dataset','donor']);st.to_csv(D/'patient_subtype_profiles.tsv',sep='\t',index=False)
mr=pd.DataFrame(matches);mr.to_csv(D/'matched_lesion_effects.tsv',sep='\t',index=False);mp=mr[mr.eligible].groupby(['branch','constraint','minpairs','dataset','donor'],as_index=False).effect.mean();mp.to_csv(D/'matched_patient_effects.tsv',sep='\t',index=False)
pd.DataFrame(detection).to_csv(D/'raw_detection_lesion.tsv',sep='\t',index=False)
def fit(v,extra=[],study=True,cov='nonrobust'):
 row=dict(n=len(v),slope=np.nan,beta=np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,status='not_estimable')
 if len(v)<4 or v.CGNL1.std()<1e-12 or v.junction.std()<1e-12:return row|dict(reason='small_n_or_no_variation')
 F=v[['CGNL1']].copy()
 if study:F=pd.concat([F,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
 for x in extra:
  if v[x].std()>1e-12:F[x]=v[x]
 F=sm.add_constant(F,has_constant='add')
 if np.linalg.matrix_rank(F.to_numpy())<F.shape[1] or len(v)-F.shape[1]<3:return row|dict(reason='rank_or_residual_df')
 m=sm.OLS(v.junction,F).fit(cov_type=cov,use_t=True);ci=m.conf_int().loc['CGNL1'];return row|dict(slope=float(m.params.CGNL1),beta=float(m.params.CGNL1*v.CGNL1.std()/v.junction.std()),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(m.pvalues.CGNL1),status='estimated',reason='',df=float(m.df_resid),covariates='|'.join(F.columns))
def save_stats(rows,name):
 t=pd.DataFrame(rows);t['FDR']=np.nan
 for fam,ix in t.groupby('family').groups.items():
  valid=t.loc[ix].index[t.loc[ix,'P'].notna()]
  if len(valid):t.loc[valid,'FDR']=multipletests(t.loc[valid,'P'],method='fdr_bh')[1]
 t['nominal_P005']=t.P<.05;t.to_csv(D/name,sep='\t',index=False);return t
models={'study':[],'study_sinusoidal':['sinusoidal_fraction'],'study_lymphatic':['lymphatic_fraction'],'study_both_fractions':['sinusoidal_fraction','lymphatic_fraction'],'study_both_depth':['sinusoidal_fraction','lymphatic_fraction','mean_logUMI']}
rows=[];rob=[]
for branch,g in pt.groupby('branch'):
 for model,extra in models.items():rows.append(dict(branch=branch,model=model,family='A04_composition_'+branch,**fit(g,extra)))
 rob.append(dict(branch=branch,model='study_both_fractions_HC3',family='A04_composition_robustness',**fit(g,models['study_both_fractions'],cov='HC3')))
cm=save_stats(rows,'composition_adjusted_models.tsv');save_stats(rob,'composition_robustness.tsv')
rows=[]
for branch in ['all','prior_reference_supported']:
 for n in ([3,5,10] if branch=='all' else [5]):
  for label in types:
   for scope in ['pooled','GSE245552','GSE285991']:
    g=st[(st.branch==branch)&(st.mincells==n)&(st.subtype==label)]
    if scope!='pooled':g=g[g.dataset==scope]
    for model,extra in [('study',[]),('study_depth',['mean_logUMI'])]:rows.append(dict(branch=branch,mincells=n,subtype=label,scope=scope,model=model,family='A04_subtype_'+branch+'_n'+str(n),**fit(g,extra)))
smr=save_stats(rows,'within_subtype_models.tsv')
def effect(v):
 x=v.effect.to_numpy();n=len(x);row=dict(n=n,positive=int((x>0).sum()),effect=float(x.mean()) if n else np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,signflip_P=np.nan)
 if n>=4 and x.std(ddof=1)>1e-12:
  ci=stats.t.interval(.95,n-1,loc=x.mean(),scale=stats.sem(x));row.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=ci[0],CI_high=ci[1])
  if n<=18:row['signflip_P']=sum(abs(np.dot(s,x)/n)>=abs(x.mean())-1e-12 for s in itertools.product([-1,1],repeat=n))/2**n
 return row
rows=[]
for branch in ['all','prior_reference_supported']:
 for constraint in ['depth_only','same_coarse_subtype']:
  for n in [3,10]:
   for scope in ['pooled','GSE245552','GSE285991']:
    g=mp[(mp.branch==branch)&(mp.constraint==constraint)&(mp.minpairs==n)]
    if scope!='pooled':g=g[g.dataset==scope]
    rows.append(dict(branch=branch,constraint=constraint,minpairs=n,scope=scope,family='A06_same_subtype_'+branch,**effect(g)))
mm=save_stats(rows,'same_subtype_matched_models.tsv');commonrows=[]
for (branch,n),g in mp.groupby(['branch','minpairs']):
 a=g[g.constraint=='depth_only'];b=g[g.constraint=='same_coarse_subtype'];common=set(a.donor)&set(b.donor)
 for label,v in [('depth_all',a),('depth_common',a[a.donor.isin(common)]),('same_subtype_common',b[b.donor.isin(common)])]:commonrows.append(dict(branch=branch,minpairs=n,version=label,n=len(v),positive=int((v.effect>0).sum()),effect=float(v.effect.mean()) if len(v) else np.nan,patients='|'.join(sorted(v.donor))))
pd.DataFrame(commonrows).to_csv(D/'matching_common_patient_diagnostic.tsv',sep='\t',index=False)
rows=[]
for scope in ['pooled','GSE245552','GSE285991']:
 g=pt[pt.branch=='all']
 if scope!='pooled':g=g[g.dataset==scope]
 for donor in g.donor:rows.append(dict(scope=scope,omitted=donor,family='A05_influence_'+scope,**fit(g[g.donor!=donor])))
save_stats(rows,'continuous_leave_one_patient_out.tsv')
# Reconcile baseline profiles and matching to previous frozen retest before interpreting.
old=pd.read_csv(B/'review_consensus_retest_20261002_v1/results/patient_continuous_profiles.tsv',sep='\t');old=old[old.stage=='cluster_primary_n10'].set_index('donor');new=pt[pt.branch=='all'].set_index('donor');assert set(old.index)==set(new.index);assert np.allclose(old.loc[new.index,['CGNL1','junction','mean_logUMI']],new[['CGNL1','junction','mean_logUMI']])
oldm=pd.read_csv(B/'review_consensus_retest_20261002_v1/results/matched_patient_effects.tsv',sep='\t');oldm=oldm[oldm.stage=='matched_primary_n10_pairs3'].set_index('donor');newm=mp[(mp.branch=='all')&(mp.constraint=='depth_only')&(mp.minpairs==3)].set_index('donor');assert set(oldm.index)==set(newm.index);assert np.allclose(oldm.loc[newm.index].effect,newm.effect)
assert all(hashlib.sha256((R/'identity'/f).read_bytes()).hexdigest()==h for f,h in freeze['hashes'].items())
q={'passed':True,'primary_patient_n':len(new),'baseline_profiles_reproduced':True,'baseline_matched_effects_reproduced':True,'identity_frozen_unchanged':True,'raw_reference_types_not_mixed':True};(R/'QA.json').write_text(json.dumps(q,indent=2))
summary={'finished':datetime.datetime.now().astimezone().isoformat(),'exit_code':0,'primary_composition':cm[cm.branch=='all'].to_dict('records'),'same_subtype_matches':mm[(mm.scope=='pooled')&(mm.constraint=='same_coarse_subtype')].to_dict('records'),'nominal_within_subtype':smr[smr.nominal_P005].to_dict('records'),'scope':'two cohorts A04 diagnosis; original overall primary negatives retained'}
(R/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
