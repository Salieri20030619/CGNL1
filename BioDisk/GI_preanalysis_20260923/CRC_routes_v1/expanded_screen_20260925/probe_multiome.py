from pathlib import Path
import anndata as ad,json,numpy as np
A=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/expanded_screen_20260925');a=ad.read_h5ad(A/'multiome/CRCLM_multiome_GEX_decontaminated.h5ad',backed='r')
print(a.shape);print(a.obs.columns.tolist());print(a.var.columns.tolist());print('layers',list(a.layers));print(a.obs.head().to_string());print(a.var.head().to_string())
counts={c:a.obs[c].value_counts().to_dict() for c in a.obs.columns if a.obs[c].nunique()<60};print(json.dumps(counts,default=str,indent=2));(A/'MULTIOME_METADATA_PROBE.json').write_text(json.dumps({'shape':a.shape,'obs_columns':list(a.obs.columns),'var_columns':list(a.var.columns),'layers':list(a.layers),'metadata_value_counts':counts},default=str,indent=2))
