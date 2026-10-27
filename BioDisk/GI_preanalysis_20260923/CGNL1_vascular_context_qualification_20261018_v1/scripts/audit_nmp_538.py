from pathlib import Path
import anndata as ad,numpy as np,h5py,json,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");s=json.loads((B/"metadata/NMP534_STATE.json").read_text());a=ad.read_h5ad(s["matrix_path"],backed="r");obs=a.obs
mapping=obs.groupby(["sample","patient","timepoint","material","report"],observed=True).size().reset_index(name="cells").to_dict("records")
ec=obs[obs.cell_type=="Endothelial cells"].groupby(["patient","timepoint","material"],observed=True).size().reset_index(name="endothelial_cells").to_dict("records")
num={}
with h5py.File(s["matrix_path"],"r") as h:
 for name in ["X","raw/X"]:
  g=h[name];d=g["data"];v=dict(stored_values=len(d),minimum=None,maximum=None,nonfinite=0,negative=0,noninteger_gt_1e6=0)
  for start in range(0,len(d),4000000):
   x=np.asarray(d[start:start+4000000],dtype=np.float64);v["nonfinite"]+=int((~np.isfinite(x)).sum());v["negative"]+=int((x<0).sum());v["noninteger_gt_1e6"]+=int((abs(x-np.rint(x))>1e-6).sum());v["minimum"]=float(x.min()) if v["minimum"] is None else min(v["minimum"],float(x.min()));v["maximum"]=float(x.max()) if v["maximum"] is None else max(v["maximum"],float(x.max()))
  num[name]=v
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),mapping=mapping,endothelial_coverage=ec,numerical=num,raw_shape=a.raw.shape,raw_target_matches=[str(g) for g in a.raw.var_names if str(g)=="CGNL1"],cell_count_discrepancy={"matrix":len(obs),"paper":75075,"difference":len(obs)-75075},mapping_validated_against_author_code=False)
(B/"results/nmp_qualification_536/AUDIT538.json").write_text(json.dumps(r,indent=2),encoding="utf-8");print(json.dumps(r));a.file.close()
