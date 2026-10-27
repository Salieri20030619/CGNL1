from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
O=Path.cwd();p=pd.read_csv(O/'results/patient_effects.tsv',sep='\t');s=pd.read_csv(O/'results/section_effects.tsv',sep='\t');e=pd.read_csv(O/'results/model_eligibility.tsv',sep='\t');b=pd.read_csv(O/'results/matching_balance.tsv',sep='\t')
note={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'post-primary diagnostic, all first-stage results known','reasons':['block E1 differs from spot main result','relative comparator E2 nominal while absolute E2 not nominal','matching support needs evaluation'],'fixed_checks':['compare spot and block E1 on exact same sections and patients','display individual patient and comparator E2 effects, no new selected marker','quantify matching support and balance; distinguish low power from extrapolation','inspect coefficient stability and rank; do not replace main result with a successful subset'],'new_primary_hypothesis':False}
(O/'POST_PRIMARY_DIAGNOSTIC_SPEC.json').write_text(json.dumps(note,indent=2),encoding='utf-8')
print('PRIMARY_PATIENTS',p[(p.variant=='primary')&(p.gene=='CGNL1')].to_string(index=False))
print('MATCH_BALANCE',b.to_string(index=False))
block=s[(s.variant=='blocks250')&(s.endpoint=='E1')&(s.gene=='CGNL1')]
same=s[(s.variant=='primary')&(s.endpoint=='E1')&(s.gene=='CGNL1')&s['sample'].isin(block['sample'])]
common=block[['sample','patient','effect']].merge(same[['sample','effect']],on='sample',suffixes=('_block','_spot'));common.to_csv(O/'results/block_same_section_comparison.tsv',sep='\t',index=False)
pp=common.groupby('patient')[['effect_block','effect_spot']].mean();pp.to_csv(O/'results/block_same_patient_comparison.tsv',sep='\t');print('EXACT SAME SUPPORT',pp.to_string());print('COMMON_MEAN',pp.mean().to_dict())
bench=p[(p.variant=='primary')&(p.endpoint=='E2')].pivot(index='patient',columns='gene',values='effect_SD');bench['comparator_median']=bench[[g for g in bench if g!='CGNL1']].median(axis=1);bench['target_minus_comparator']=bench.CGNL1-bench.comparator_median;bench.to_csv(O/'results/E2_individual_specificity.tsv',sep='\t');print('E2_SPECIFICITY',bench.to_string())
print('ELIGIBILITY',e.groupby(['variant','endpoint','status']).size().to_string())
diag={'same_section_block_mean':float(pp.effect_block.mean()),'same_section_spot_mean':float(pp.effect_spot.mean()),'same_support_patients':len(pp),'balanced_matched_sections_ge20':int(((b.pairs>=20)&b.balanced).sum()),'matched_sections_ge20':int((b.pairs>=20).sum()),'matching_not_guaranteed_balanced':True,'no_new_inferential_tests':True}
(O/'DIAGNOSTIC_COMPLETE.json').write_text(json.dumps(diag,indent=2),encoding='utf-8');print(json.dumps(diag))
