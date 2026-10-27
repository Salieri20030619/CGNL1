from pathlib import Path
import json,numpy as np,pandas as pd,statsmodels.api as sm
from scipy import sparse
r=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/review_consensus_retest_20261002_v1');d=r/'results'
new=pd.read_csv(d/'patient_continuous_profiles.tsv',sep='\t');old=pd.read_csv(r.parent/'threshold_ladder_20261002_v1/results/patient_profiles_all_stages.tsv',sep='\t');old=old[old.dataset.isin(['GSE245552','GSE285991'])]
def fit(v,depth=False,cov='nonrobust'):
 F=pd.concat([v[['CGNL1']],pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1)
 if depth:F['depth']=v.mean_logUMI
 F=sm.add_constant(F,has_constant='add');m=sm.OLS(v.junction,F).fit(cov_type=cov,use_t=True)
 return dict(n=len(v),slope=m.params.CGNL1,beta=m.params.CGNL1*v.CGNL1.std()/v.junction.std(),P=m.pvalues.CGNL1,df=m.df_resid)
rows=[];prim=new[new.stage=='cluster_primary_n10']
for stage in ['B0_broad','Q3_mito25','I1_markers2','N2_cells50']:
 ov=old[old.stage==stage];assert not ov.donor.duplicated().any();common=set(ov.donor)&set(prim.donor)
 for label,v in [('old_all',ov),('old_common',ov[ov.donor.isin(common)]),('reviewed_common',prim[prim.donor.isin(common)])]:rows.append(dict(comparison_stage=stage,version=label,purpose='posthoc_identity_diagnostic',**fit(v)))
pd.DataFrame(rows).to_csv(d/'identity_same_cohort_common_patient_diagnostic.tsv',sep='\t',index=False)
a=new[new.stage=='cluster_primary_n50'];common=set(a.donor);steps=[]
for stage,v in [('n10_all',prim),('n10_fixed_n50_patients',prim[prim.donor.isin(common)]),('n50_qualifying_lesions',a)]:
 for depth in [False,True]:steps.append(dict(stage=stage,depth_adjusted=depth,**fit(v,depth)))
pd.DataFrame(steps).to_csv(d/'n50_common_patient_diagnostic.tsv',sep='\t',index=False)
loo=[dict(omitted='none_HC3',**fit(a,True,'HC3'))]+[dict(omitted=x,**fit(a[a.donor!=x],True)) for x in a.donor];pd.DataFrame(loo).to_csv(d/'n50_depth_posthoc_robustness.tsv',sep='\t',index=False)
mr=pd.read_csv(d/'matched_lesion_results.tsv',sep='\t');m0=mr[(mr.stage=='unmatched_n10')&mr.eligible];m1=mr[(mr.stage=='matched_primary_n10_pairs3')&mr.eligible];assert set(m0['sample'])==set(m1['sample'])
qa=[];detect=[];q=pd.read_csv(d/'lesion_continuous_profiles.tsv',sep='\t')
for rec in json.loads((r/'CACHE_INDEX.json').read_text()):
 y=sparse.load_npz(rec['path']);genes=pd.read_csv(rec['genes'],sep='\t').gene;g=pd.read_csv(rec['cells'],sep='\t');assert y.shape==(len(genes),len(g));counts=y[np.flatnonzero(genes.eq('CGNL1'))].toarray().ravel();lib=np.asarray(y.sum(axis=0)).ravel();pos=counts>0
 if rec['dataset']=='GSE245552':assert g.mito_QC_current.isna().all() and g.mito_retest.isna().all()
 detect.append(dict(dataset=rec['dataset'],donor=rec['donor'],sample=rec['sample'],cells=len(g),CGNL1_counts=counts.sum(),CGNL1_detected_cells=int(pos.sum()),CGNL1_detection_fraction=float(pos.mean()) if len(g) else np.nan,positive_mean_logUMI=np.log1p(lib[pos]).mean() if pos.any() else np.nan,negative_mean_logUMI=np.log1p(lib[~pos]).mean() if (~pos).any() else np.nan))
 if len(g)>=10:
  ix=[int(np.flatnonzero(genes.eq(x))[0]) for x in ['CGNL1','CLDN5','ESAM','TJP1','OCLN']];p=np.log1p(np.asarray(y[ix].sum(axis=1)).ravel()/lib.sum()*1e4);z=q[(q.stage=='cluster_primary_n10')&(q['sample']==rec['sample'])].iloc[0];assert np.allclose([p[0],p[1:].mean()],[z.CGNL1,z.junction]);qa.append(rec['sample'])
pd.DataFrame(detect).to_csv(d/'CGNL1_detection_and_depth_diagnostic.tsv',sep='\t',index=False)
t=pd.read_csv(d/'continuous_models.tsv',sep='\t');assert len(t)==108;assert np.isclose(fit(prim)['P'],t[t.is_primary&(t.scope=='pooled')].P.iloc[0])
Q={'status':'passed','raw_cache_libraries':27,'reconstructed_primary_lesions':len(qa),'continuous_model_rows':len(t),'matched_rows':len(pd.read_csv(d/'matched_models.tsv',sep='\t')),'unmatched_and_matched_same_eligible_lesions':sorted(m0['sample'].tolist()),'D02_mito_missing_preserved':True,'posthoc_diagnostics_do_not_replace_primary':True};(r/'QA.json').write_text(json.dumps(Q,indent=2))
for fn in ['identity_same_cohort_common_patient_diagnostic.tsv','n50_common_patient_diagnostic.tsv','n50_depth_posthoc_robustness.tsv']:
 print(fn);print(pd.read_csv(d/fn,sep='\t').to_string(index=False))
print(json.dumps(Q))
# FDR for the three posthoc diagnostic tables was recorded separately on completion.

from statsmodels.stats.multitest import multipletests
for f in ['identity_same_cohort_common_patient_diagnostic.tsv','n50_common_patient_diagnostic.tsv','n50_depth_posthoc_robustness.tsv']:
 t=pd.read_csv(d/f,sep='\t');t['family']='posthoc_'+f.removesuffix('.tsv');t['FDR']=multipletests(t.P,method='fdr_bh')[1];t['interpretation']='diagnostic_only_not_independent_confirmation';t.to_csv(d/f,sep='\t',index=False)
