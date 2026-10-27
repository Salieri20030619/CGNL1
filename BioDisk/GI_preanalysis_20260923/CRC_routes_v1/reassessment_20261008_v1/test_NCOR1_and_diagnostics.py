from pathlib import Path
import json,datetime,hashlib,gzip
import pandas as pd,numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1');D=R/'perturbation_inputs';O=R/'results';spec=json.loads((R/'REASSESSMENT_SPEC.json').read_text())
p=D/'GSE237653_2023_Teichmann_NCoR1_RNAseq_NormCounts.xlsx';t=pd.read_excel(p,sheet_name='NormCounts',engine='openpyxl');cols=['CTL_1','CTL_2','CTL_3','siNCoR_1','siNCoR_2','siNCoR_3'];assert set(cols)<=set(t.columns);assert t['Ensembl gene id'].is_unique
genes=['CGNL1','NCOR1','DLL4','NOTCH4']+spec['R4_real_perturbation']['positive_controls'];rows=[];vals=[]
for gene in genes:
 mask=t['Ensembl gene'].eq(gene)
 if gene=='CGNL1':mask |= t['Ensembl gene id'].eq('ENSG00000128849')
 hit=t[mask]
 if len(hit)!=1:
  rows.append(dict(gene=gene,status='missing_or_ambiguous_mapping',hits=len(hit)));continue
 v=hit.iloc[0][cols].to_numpy(float);assert np.isfinite(v).all() and (v>=0).all();a,b=v[:3],v[3:]
 rows.append(dict(dataset='GSE237653',cell_type='HUVEC',gene=gene,ensembl=hit.iloc[0]['Ensembl gene id'],status='descriptive_normalized_counts',control_mean=float(a.mean()),treated_mean=float(b.mean()),log2_ratio_means_plus1=float(np.log2((b.mean()+1)/(a.mean()+1))),treated_above_control_mean=int((b>a.mean()).sum()),all_treated_above_all_controls=bool(b.min()>a.max()),all_treated_below_all_controls=bool(b.max()<a.min()),P=np.nan,FDR=np.nan,reason_no_new_P='normalized_counts_only; per-column donor/lot pairing not established; avoid new donor-level inference'))
 for col,val in zip(cols,v):vals.append(dict(gene=gene,column=col,condition='control' if col.startswith('CTL') else 'siNCOR1',normalized_count=float(val)))
pd.DataFrame(rows).to_csv(O/'R4_NCOR1_target_and_controls.tsv',sep='\t',index=False);pd.DataFrame(vals).to_csv(O/'R4_NCOR1_individual_values.tsv',sep='\t',index=False)
(O/'R4_NCOR1_QUALIFICATION.json').write_text(json.dumps({'input':'author NormCounts, not raw count input despite integer-looking values','samples':cols,'metadata':'HUVEC lots405Z013,408Z014,416Z042 listed collectively; exact column-to-lot mapping not specified','treatment':'siNCOR1 versus negative siRNA;48h','new_P_values':False,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
# Influence diagnostics on every estimable R2 model, not only nominally positive branches.
pt=pd.read_csv(O/'R2_patient_lineage_profiles.tsv',sep='\t');models=pd.read_csv(O/'R2_state_models.tsv',sep='\t');diag=[];lo=[];verify=[]
for _,r in models[models.variance.eq('OLS')].iterrows():
 v=pt[pt.lineage.eq(r.lineage)];v=v if r.scope=='pooled' else v[v.dataset.eq(r.scope)];v=v.dropna(subset=['CGNL1','mean_logUMI','hypoxia',r.program]).reset_index(drop=True)
 Z=pd.DataFrame({'const':np.ones(len(v)),'depth':v.mean_logUMI,'hypoxia':v.hypoxia});Z=pd.concat([Z,pd.get_dummies(v.dataset,drop_first=True,dtype=float)],axis=1);F=Z.copy();F['CGNL1']=v.CGNL1;m=sm.OLS(v[r.program],F).fit();inf=m.get_influence();rx=sm.OLS(v.CGNL1,Z).fit();VIF=1/(1-rx.rsquared) if rx.rsquared<1 else np.inf
 diag.append(dict(scope=r.scope,lineage=r.lineage,program=r.program,n=len(v),target_VIF=VIF,max_leverage=float(inf.hat_matrix_diag.max()),max_Cooks=float(inf.cooks_distance[0].max()),max_Cooks_patient=v.patient.iloc[np.argmax(inf.cooks_distance[0])]))
 verify.append(dict(scope=r.scope,lineage=r.lineage,program=r.program,old_P=r.P,recomputed_P=float(m.pvalues.CGNL1),passed=bool(np.isclose(m.pvalues.CGNL1,r.P))))
 for i in range(len(v)):
  keep=np.arange(len(v))!=i;mm=sm.OLS(v.loc[keep,r.program],F.loc[keep]).fit();lo.append(dict(scope=r.scope,lineage=r.lineage,program=r.program,omitted_patient=v.patient.iloc[i],slope=float(mm.params.CGNL1),P=float(mm.pvalues.CGNL1),role='influence_diagnostic_not_separate_discovery'))
pd.DataFrame(diag).to_csv(O/'R2_design_influence_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(lo).to_csv(O/'R2_leave_one_patient_influence.tsv',sep='\t',index=False);pd.DataFrame(verify).to_csv(O/'R2_model_numeric_verification.tsv',sep='\t',index=False);assert all(r['passed'] for r in verify)
fdrchecks=[]
for filename in ['R2_state_models.tsv','R3_paired_summary.tsv']:
 t=pd.read_csv(O/filename,sep='\t')
 for fam,g in t.groupby('family'):
  valid=g.P.notna();expected=multipletests(g.loc[valid,'P'],method='fdr_bh')[1];ok=np.allclose(expected,g.loc[valid,'FDR']);fdrchecks.append(dict(file=filename,family=fam,n=int(valid.sum()),passed=bool(ok)))
assert all(r['passed'] for r in fdrchecks);pd.DataFrame(fdrchecks).to_csv(O/'FDR_numeric_verification.tsv',sep='\t',index=False)
(O/'NCOR1_DIAGNOSTICS_COMPLETE.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'NCOR1_descriptive_completed':True,'R2_models_verified':len(verify),'FDR_families_verified':len(fdrchecks),'new_donor_inference_for_cell_culture':False},indent=2));print(pd.DataFrame(rows).query("gene=='CGNL1' or gene=='NCOR1'").to_json(orient='records'))
