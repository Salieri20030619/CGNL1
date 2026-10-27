from pathlib import Path
import pandas as pd,fitz
p=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation\final_exploration_20260924_v2\F01_SOURCE_STATE\donor_state_effects.tsv');t=pd.read_csv(p,sep='\t');t=t[(t.target=='CGNL1')&(t.lineage=='endothelial')];print(t.to_string(index=False))
p=next(Path('D:/').glob('电脑管家迁移文件/xwechat_files/wxid_913cjq0lvmfd32_3e57/msg/file/2026-09/*v1.1.pdf'));d=fitz.open(p);d[1].get_pixmap().save(r'C:\Users\31937\CGNL1_review_consensus_retest_20261002\plan_pdf_audit\page_2.png')