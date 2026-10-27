from pathlib import Path
import fitz
p=Path(r'D:\电脑管家迁移文件\xwechat_files\wxid_913cjq0lvmfd32_3e57\msg\file\2026-09\CGNL1纯生信研究计划_客户版_v1.1.pdf');d=fitz.open(p);out=Path(r'C:\Users\31937\CGNL1_review_consensus_retest_20261002\plan_pdf_audit');out.mkdir(exist_ok=True)
texts=[p.get_text() for p in d];(out/'extracted.txt').write_text('\n'.join('PAGE '+str(i+1)+'\n'+t for i,t in enumerate(texts)),encoding='utf-8');print('PAGES',len(d))
for i,t in enumerate(texts):
 if any(s in t for s in ['A03','A04','A05','A06','已有证据决定','外部内皮状态','执行顺序']):
  print('\nPAGE',i+1,'\n',t);d[i].get_pixmap(matrix=fitz.Matrix(1,1)).save(str(out/('page_'+str(i+1)+'.png')))