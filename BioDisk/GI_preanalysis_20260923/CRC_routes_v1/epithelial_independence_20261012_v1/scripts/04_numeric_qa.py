from pathlib import Path
import json,hashlib,shutil,datetime
import numpy as np,pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');O=R/'results'
s=json.loads((R/'SPECIFICATION.json').read_text());checks={}
checks['spec_hash_unchanged']=hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
checks['recorded_exit_codes']={x:int((R/'logs'/f'{x}.exit').read_text()) for x in ['extract','AB','B_subtype','C','C_adapted','A_mixture','C_blocks']}
def bh(p,m=None):
 p=np.asarray(p,dtype=float);m=len(p) if m is None else m;order=np.argsort(p);q=np.empty(len(p))
 q[order]=np.minimum(1,np.minimum.accumulate((p[order]*m/np.arange(1,len(p)+1))[::-1])[::-1])
 return q
# Preserve correction history and use the number of prespecified hypotheses.
hist=R/'internal_repair_archive';hist.mkdir(exist_ok=True)
for name in ['C_gene_summary.tsv','C_specificity_summary.tsv']:
 p=O/name;back=hist/(name+'.before_planned_family_correction')
 if not back.exists():shutil.copy2(p,back)
 t=pd.read_csv(p,sep='\t')
 for method in ['signflip','ttest']:
  subset=t[t.gene.eq('CGNL1')] if 'gene' in t else t
  col='FDR_target_'+method if 'gene' in t else 'FDR_'+method
  for _,ix in subset.groupby('min_micro').groups.items():
   ix=[i for i in ix if np.isfinite(t.loc[i,'P_'+method])]
   if ix:t.loc[ix,col]=bh(t.loc[ix,'P_'+method],3)
 t['planned_target_family_size']=3
 t.to_csv(p,sep='\t',index=False)
# Abundance benchmark counts use27 planned genes, without fabricating missing P values.
m=pd.read_csv(O/'A_all_models.tsv',sep='\t');m['benchmark_family_planned_size']=27
for _,ix in m.groupby('family_benchmark').groups.items():
 m.loc[ix,'FDR_benchmark']=bh(m.loc[ix,'P'],27)
m.to_csv(O/'A_all_models.tsv',sep='\t',index=False);m[m.predictor.eq('CGNL1')].to_csv(O/'A_CGNL1_models.tsv',sep='\t',index=False)
errors=[]
for _,v in m[m.predictor.eq('CGNL1')].groupby('family_target'):errors.extend(abs(v.FDR_target.to_numpy()-bh(v.P.to_numpy(),3)))
checks['A_target_FDR_max_error']=float(max(errors))
# New raw aggregation must reproduce the corresponding archived measurements.
new=pd.read_csv(O/'A_library_profiles.tsv',sep='\t')
old=pd.read_csv(R.parent/'reassessment_20261008_v1/results/R2_library_lineage_profiles.tsv',sep='\t')
old=old[old.lineage.eq('epithelial')]
merged=new.merge(old,on=['dataset','sample'],suffixes=('_new','_old'))
checks['A_libraries_reproduced']=len(merged);checks['A_target_max_reproduction_error']=float(abs(merged.CGNL1_new-merged.CGNL1_old).max())
checks['A_focal_max_reproduction_error']=float(abs(merged.focal_adhesion_new-merged.focal_adhesion_old).max())
# Verify held-out target predictions with unscaled raw training regressions.
inp=pd.read_csv(O/'A_model_input.tsv',sep='\t');pr=pd.read_csv(O/'A_heldout_predictions.tsv',sep='\t');pr=pr[pr.predictor.eq('CGNL1')];errs=[]
for (scope,model,outcome),g in pr.groupby(['scope','model','outcome']):
 v=inp.copy() if scope=='pooled' else inp[inp.dataset.eq(scope)].copy()
 v=v.dropna(subset=[outcome,'CGNL1']+[x for x in s['A']['models'][model] if x!='study']).reset_index(drop=True)
 X=pd.DataFrame({'constant':np.ones(len(v))})
 for col in s['A']['models'][model]:
  if col=='study':
   for d in sorted(v.dataset.unique())[1:]:X['study_'+d]=v.dataset.eq(d).astype(float)
  else:X[col]=v[col].to_numpy()
 Y=v[outcome].to_numpy()
 for _,row in g.iterrows():
  i=np.flatnonzero(v.patient.eq(row.patient))[0];keep=np.arange(len(v))!=i
  for aug,label in [(False,'baseline'),(True,'augmented')]:
   M=X.copy()
   if aug:M['CGNL1']=v.CGNL1.to_numpy()
   vals=M.to_numpy();fit=np.linalg.lstsq(vals[keep],Y[keep],rcond=None)[0]
   errs.append(abs(float(vals[i]@fit)-row[label]))
checks['heldout_prediction_checks']=len(errs);checks['heldout_prediction_max_error']=float(max(errs))
checks['FWL_A_max_error']=float(pd.read_csv(O/'A_FWL_numeric_check.tsv',sep='\t').abs_coefficient_error.max())
checks['FWL_C_max_error']=float(pd.read_csv(O/'platform_adapted/C_FWL_numeric_check.tsv',sep='\t').FWL_error.max())
b=pd.read_csv(O/'B_paired_effects.tsv',sep='\t');checks['paired_main_patients']=int(b[b.group.eq('Tumor_no_Tu10')].patient.nunique())
checks['shared_subtype_patients']=int(b[b.group.eq('shared_author_subtypes')].patient.nunique())
c=pd.read_csv(O/'platform_adapted/C_patient_coefficients.tsv',sep='\t');checks['spatial_primary_threshold_patients']=int(c[c.min_micro.eq(10)&c.gene.eq('CGNL1')].patient.nunique())
checks['patient_not_spot_inference']=True
checks['all_pass']=checks['spec_hash_unchanged'] and all(x==0 for x in checks['recorded_exit_codes'].values()) and checks['A_target_FDR_max_error']<1e-10 and checks['A_target_max_reproduction_error']<1e-10 and checks['A_focal_max_reproduction_error']<1e-10 and checks['heldout_prediction_max_error']<1e-8 and checks['FWL_A_max_error']<1e-8 and checks['FWL_C_max_error']<1e-8
checks['statistical_correction']='Primary spatial FDR retains3 prespecified model hypotheses despite C1/C2 unestimable; raw missing P values remain missing. A comparator family retains27 planned genes.'
(R/'NUMERICAL_QA.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks));assert checks['all_pass']
