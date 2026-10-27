from pathlib import Path
p=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/scripts/figures_verify.py')
s=p.read_text();s=s.replace('import matplotlib.pyplot as plt','import matplotlib.pyplot as plt\nfrom scipy.stats import t as student_t');s=s.replace('xerr=z.SE*1.96','xerr=z.SE*student_t.ppf(.975,15)');s=s.replace('xerr=z.SE_composition*1.96','xerr=z.SE_composition*student_t.ppf(.975,13)');s=s.replace('approximate 95% CI','95% t CI');s=s.replace('diagnostic gene panels','exploratory gene panels');p.write_text(s)
