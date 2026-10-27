from pathlib import Path
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');p=R/'scripts'/'client_figures.py';s=p.read_text().replace("'font.size':11","'font.size':14").replace('fontsize=12','fontsize=16').replace('fontsize=13','fontsize=17').replace('fontsize=10','fontsize=13').replace('fontsize=9','fontsize=12').replace('fontsize=8','fontsize=11');p.write_text(s)
