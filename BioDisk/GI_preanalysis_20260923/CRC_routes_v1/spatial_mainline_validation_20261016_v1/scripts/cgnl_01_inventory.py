from pathlib import Path
import json,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');O=Path.cwd()
folders=['reassessment_20261008_v1','feasibility_audit_20261008_v1','epithelial_independence_20261012_v1','spatial_v3']
inventory=[]
for f in folders:
 for p in (R/f).rglob('*'):
  if p.is_file() and p.suffix.lower() in ['.py','.r','.json','.parquet','.h5ad','.npz','.pkl']:
   inventory.append({'path':str(p),'bytes':p.stat().st_size})
(O/'metadata/EXISTING_INPUT_INVENTORY.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
for x in inventory:
 if x['path'].endswith(('.py','.r','.parquet','.h5ad','.npz','.pkl')) and ('spatial_v3' not in x['path'] or x['path'].endswith(('.parquet','.npz','.pkl'))):print(x)
for rel in ['epithelial_independence_20261012_v1/scripts/03_analyze_C.py','epithelial_independence_20261012_v1/scripts/03b_spatial_platform_adapted.py']:
 print('\nSOURCE',rel);print((R/rel).read_text(encoding='utf-8'))
