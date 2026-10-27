from pathlib import Path
import json,re
O=Path.cwd()/'mechanism_investigation_v2';text=(O/'literature/PMGC2026_accepted.txt').read_text(encoding='utf-8');pages=text.split('\f');terms=['spatial transcript','mesothelial','mesotheli','quality control','mitochondrial','infercnv','copykat','10x','visium'];hits=[]
for i,p in enumerate(pages):
 if any(t in p.lower() for t in terms):hits.append({'page':i+1,'text':p})
(O/'metadata/PMGC2026_DESIGN_143.json').write_text(json.dumps({'pages':hits,'total_pdf_text_pages':len(pages)},ensure_ascii=False,indent=2),encoding='utf-8')
for h in hits:
 if h['page']>=20 and h['page']<=30:print('PAGE',h['page'],h['text'])
