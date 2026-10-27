from pathlib import Path
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');p=R/'scripts'/'client_figures.py';s=p.read_text();s=s.replace('10⁻⁷','10^-7').replace('10⁻⁹','10^-9').replace('10⁻⁸','10^-8').replace('log₂','log2');s=s.replace("'font.size':11","'font.size':14");p.write_text(s)
