from pathlib import Path
import anndata as ad,numpy as np,json,datetime
from scipy import sparse
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1")
spec=dict(scope="Post hoc marker-only phenotype sensitivity across ALL nuclei, not restricted to author endothelial cluster. Not a validated cell classifier or HGP test.",pan_markers=["PECAM1","CDH5","CLDN5","ESAM"],pan_minimum=2,lymphatic_markers=["PROX1","PDPN"],positive_tolerance=1e-6,lymphatic_rule="pan minimum AND both lymphatic markers positive",comparator_rule="pan minimum AND both lymphatic markers absent at numerical tolerance; lymphatic-marker-negative comparator, NOT definitive blood endothelial identity",excluded_rule="one lymphatic marker positive is ambiguous",normalization="Per-nucleus CP10k denominator excludes CGNL1",coverage_thresholds=[1,3,5,10],biological_P=None)
(B/"metadata/MARKER_SENSITIVITY_SPEC_532.json").write_text(json.dumps(spec,indent=2),encoding="utf-8")
state=json.loads((B/"metadata/PROCESSED525_STATE.json").read_text());a=ad.read_h5ad(state["verified_path"]);m=sparse.csr_matrix(a.X,dtype=np.float64)
def gene(g):return m[:,a.var_names.get_loc(g)].toarray().ravel()
pan=sum((gene(g)>1e-6).astype(int) for g in spec["pan_markers"])>=2
l1=gene("PROX1")>1e-6;l2=gene("PDPN")>1e-6
masks={"lymphatic_marker_positive":pan&l1&l2,"lymphatic_marker_negative":pan&~l1&~l2,"ambiguous":pan&(l1^l2)}
t=gene("CGNL1");denom=np.asarray(m.sum(axis=1)).ravel()-t;norm=np.divide(t*10000,denom,out=np.zeros(len(t)),where=denom>0);patient=np.asarray(a.obs["Patient"].astype(str));ctype=np.asarray(a.obs["Cell_type"].astype(str));rows=[]
for p in sorted(set(patient)):
 for name,mask in masks.items():
  s=mask&(patient==p);rows.append(dict(patient=p,phenotype=name,nuclei=int(s.sum()),mean_CP10k=float(norm[s].mean()) if s.any() else None,positive=int((t[s]>1e-6).sum()),author_cell_types={str(x):int(y) for x,y in zip(*np.unique(ctype[s],return_counts=True))}))
lookup={(x["patient"],x["phenotype"]):x for x in rows};contrasts=[]
for minimum in spec["coverage_thresholds"]:
 pairs=[]
 for p in sorted(set(patient)):
  x=lookup[p,"lymphatic_marker_positive"];y=lookup[p,"lymphatic_marker_negative"]
  if min(x["nuclei"],y["nuclei"])>=minimum:pairs.append(dict(patient=p,n_a=x["nuclei"],n_b=y["nuclei"],delta_CP10k=x["mean_CP10k"]-y["mean_CP10k"],delta_detection_fraction=x["positive"]/x["nuclei"]-y["positive"]/y["nuclei"]))
 contrasts.append(dict(minimum=minimum,patients=len(pairs),pairs=pairs))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),spec=spec,rows=rows,contrasts=contrasts,mainline_validated=False)
out=B/"results/independent_qualification_526/MARKER_SENSITIVITY532.json";out.write_text(json.dumps(result,indent=2),encoding="utf-8");print(json.dumps(dict(totals={k:int(v.sum()) for k,v in masks.items()},contrasts=contrasts)))
