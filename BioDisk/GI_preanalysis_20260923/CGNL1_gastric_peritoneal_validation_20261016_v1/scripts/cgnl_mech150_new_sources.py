from pathlib import Path
import json,zipfile,hashlib,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';p=O/'report/review_update_v2/spatial_results_bundle.zip';dest=O/'results/PMGC2026_spatial';dest.mkdir(exist_ok=True)
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None;m=json.loads(z.read('MANIFEST.json'))
 for n,h in m.items():
  b=z.read(n);assert hashlib.sha256(b).hexdigest()==h;q=(dest/n).resolve();assert q.is_relative_to(dest.resolve())
  if q.exists():assert q.read_bytes()==b
  else:q.write_bytes(b)
print('Spatial output source hashes verified')
a=pd.read_csv(O/'results/GSE308231_source_profiles.tsv',sep='\t');q=pd.read_csv(O/'results/GSE308231_QC.tsv',sep='\t');print('QC',q.to_string(index=False))
for scheme in ['primary','sensitivity']:
 v=a[a.scheme==scheme];print('SCHEME',scheme);print(v[['sample','type','group','n_cells','target_UMI','target_detection','CGNL1_ln1pCP10k']].to_string(index=False))
print('QA119', (O/'results/GSE308231_SOURCE_QA_119.json').exists())
