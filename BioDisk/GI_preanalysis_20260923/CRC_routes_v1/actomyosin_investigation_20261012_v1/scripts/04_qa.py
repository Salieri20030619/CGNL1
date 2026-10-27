from common import *
import importlib.util
checks={'spec_unchanged':hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']}
checks['exit_codes']={x:int((R/'logs'/f'{x}.exit').read_text()) for x in ['components','external','knockout']}
t=pd.read_csv(O/'D1_D2_component_models.tsv',sep='\t')
prior=pd.read_csv(R.parent/'epithelial_independence_20261012_v1/results/A_CGNL1_models.tsv',sep='\t')
old=prior[prior.scope.eq('pooled')&prior.model.eq('M1')&prior.outcome.eq('actomyosin')&prior.variance.eq('HC3')].iloc[0]
new=t[t.scope.eq('pooled')&t.model.eq('M1')&t.outcome.eq('full')].iloc[0]
checks['prior_effect_reproduced']=abs(new.beta-old.beta_standardized)<1e-10 and abs(new.P-old.P)<1e-10
v=pd.read_csv(O/'D1_D3_model_input.tsv',sep='\t');full,p,X=fitted(v,S['models']['M1'],'score_full');decomp=pd.read_csv(O/'D1_effect_decomposition.tsv',sep='\t')
checks['additive_coefficient_error']=float(abs(decomp.contribution_to_full_slope.sum()-full['slope']))
err=[]
for _,ix in t.groupby('family').groups.items():
 fam=t.loc[ix[0],'family'];n=7 if fam.endswith('7') else 3 if fam.endswith('3') else 1
 err.extend(abs(t.loc[ix,'FDR']-bh(t.loc[ix,'P'],n)))
checks['component_FDR_error']=float(max(err))
ext=pd.read_csv(O/'D4_patient_profiles.tsv',sep='\t');checks['external_author_patients']=int(ext[ext.definition.eq('author_epithelial_QC')].patient.nunique())
checks['external_score_reconstruction_error']=float(abs(ext[S['genes']].mean(1)-ext.actomyosin).max())
# Clustered endpoint degrees of freedom audited explicitly.
cluster=json.loads((O/'D4_within_subtype.json').read_text());checks['subtype_patient_clusters']=cluster['patients'];checks['subtype_inference_df']=cluster.get('effective_inference_df')
ko=pd.read_csv(O/'D5_KO_individual_gene_values.tsv',sep='\t')
checks['KO_sample_count']=int(ko['sample'].nunique());checks['KO_expression_nonnegative']=bool((ko.CPM>=0).all());checks['KO_no_fabricated_P']=not any(x in ko.columns for x in ['P','pvalue'])
qual=json.loads((O/'D5_KO_QUALIFICATION.json').read_text());checks['KO_gene_coverage']=len(qual['core_symbols_used']);checks['KO_single_clone_disclosed']=qual['KO_clones_per_genotype']==1
checks['all_pass']=checks['spec_unchanged'] and all(x==0 for x in checks['exit_codes'].values()) and checks['prior_effect_reproduced'] and checks['additive_coefficient_error']<1e-8 and checks['component_FDR_error']<1e-10 and checks['external_score_reconstruction_error']<1e-10 and checks['KO_sample_count']==12 and checks['KO_expression_nonnegative'] and checks['KO_no_fabricated_P'] and checks['KO_single_clone_disclosed']
checks={k:(v.item() if isinstance(v,np.generic) else v) for k,v in checks.items()}
(R/'NUMERICAL_QA.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks));assert checks['all_pass']
