from pathlib import Path
import anndata as ad,h5py,json,numpy as np,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1")
s=json.loads((B/"metadata/NMP534_STATE.json").read_text());assert s["stage"]=="verified_single_cell_ready_for_qualification" and s["zip_md5_pass"] and s["single_cell_zip_crc_pass"]
a=ad.read_h5ad(s["matrix_path"],backed="r")
obs={}
for col in a.obs.columns:
 values=a.obs[col].astype(str);n=int(values.nunique());obs[col]={"unique_count":n,"value_counts":{str(k):int(v) for k,v in values.value_counts().items()} if n<=60 else None,"example_values":values.head(3).tolist()}
features={"var_columns":list(a.var.columns),"target_var_names":[str(x) for x in a.var_names if str(x).upper()=="CGNL1"],"target_other_fields":{c:a.var.index[a.var[c].astype(str).str.upper()=="CGNL1"].astype(str).tolist() for c in a.var.columns}}
with h5py.File(s["matrix_path"],"r") as h:
 def info(g):
  if isinstance(g,h5py.Dataset):return dict(shape=g.shape,dtype=str(g.dtype))
  return dict(encoding=str(g.attrs.get("encoding-type")),children={k:dict(shape=g[k].shape,dtype=str(g[k].dtype)) for k in g if isinstance(g[k],h5py.Dataset)})
 storage={"X":info(h["X"]),"layers":{k:info(h["layers"][k]) for k in h.get("layers",{})},"raw_present":"raw" in h and isinstance(h["raw"],h5py.Group)}
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),shape=a.shape,obs=obs,features=features,storage=storage,uns_keys=list(a.uns),source=s["matrix_path"],measurement_scale_qualified=False,mainline_validated=False)
p=B/"results/nmp_qualification_536";p.mkdir(exist_ok=True);(p/"METADATA536.json").write_text(json.dumps(r,indent=2),encoding="utf-8");a.file.close();print(json.dumps(r))
