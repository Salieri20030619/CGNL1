from pathlib import Path
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');s=(R/'scripts'/'state_correlations.py').read_text().replace("'independent_scRNA'","'independent_scRNA_refined'").replace("'state_correlation_plan.json'","'state_correlation_plan_refined.json'");(R/'scripts'/'state_correlations_refined.py').write_text(s)
