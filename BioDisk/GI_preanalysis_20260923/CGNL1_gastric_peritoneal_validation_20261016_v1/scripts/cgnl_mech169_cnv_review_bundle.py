from pathlib import Path
import json,zipfile,hashlib,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';root=O/'results/CNV_identity';complete=json.loads((O/'results/CNV_IDENTITY_COMPLETE_167.json').read_text(encoding='utf-8'));assert complete['all_pass'];a=pd.read_csv(root/'group_descriptors.tsv',sep='\t');b=pd.read_csv(root/'mesothelial_profile_similarity.tsv',sep='\t');print(a[a.group.isin(['immune_like','endothelial_like','strict_mesothelial','conservative_epithelial'])].to_string(index=False));print('PROFILE_CORRELATIONS',b.to_string(index=False));files=list(root.glob('*/*heatmap.png'))+[root/'group_descriptors.tsv',root/'mesothelial_profile_similarity.tsv']+list(root.glob('*/CASE_COMPLETE.json'));dest=O/'report/CNV_review_bundle.zip';manifest={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,str(p.relative_to(root)))
 z.writestr('MANIFEST.json',json.dumps(manifest,indent=2).encode('utf-8'))
print('BUNDLE',str(dest),dest.stat().st_size,hashlib.sha256(dest.read_bytes()).hexdigest())
