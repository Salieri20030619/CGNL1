from pathlib import Path
import anndata as ad,numpy as np,json,itertools,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");s=json.loads((B/"metadata/NMP534_STATE.json").read_text());a=ad.read_h5ad(s["matrix_path"],backed="r");mask=(a.obs.cell_type=="Endothelial cells")&(a.obs.timepoint=="T0");c=a[mask].to_memory();a.file.close();m=c.raw.X.astype(np.float64);target=c.raw.var_names.get_loc("CGNL1");t=m[:,target].toarray().ravel();rows=[]
for patient in sorted(set(c.obs.patient)):
 for tissue in ["CRLM","adj.liver"]:
  ix=np.asarray((c.obs.patient==patient)&(c.obs.material==tissue));total=float(m[ix].sum());v=float(t[ix].sum());denom=total-v;rows.append(dict(patient=patient,tissue=tissue,n_cells=int(ix.sum()),CGNL1_count=v,positive_cells=int((t[ix]>0).sum()),non_target_UMI=denom,log2CPM1=float(np.log2(1+1e6*v/denom)) if denom>0 else None))
contrasts=[]
for cutoff in [20,10,50]:
 pairs=[]
 for patient in sorted(set(c.obs.patient)):
  x=next(r for r in rows if r["patient"]==patient and r["tissue"]=="CRLM");y=next(r for r in rows if r["patient"]==patient and r["tissue"]=="adj.liver")
  if min(x["n_cells"],y["n_cells"])>=cutoff:pairs.append(dict(patient=patient,delta=x["log2CPM1"]-y["log2CPM1"]))
 d=np.array([r["delta"] for r in pairs]);valid=len(d)>1 and np.any(d!=0);pvalue=None
 if valid:
  perm=np.array([np.mean(d*np.array(sign)) for sign in itertools.product([-1,1],repeat=len(d))]);pvalue=float(np.mean(abs(perm)>=abs(d.mean())-1e-12))
 contrasts.append(dict(cutoff=cutoff,n=len(d),pairs=pairs,mean=float(d.mean()) if len(d) else None,median=float(np.median(d)) if len(d) else None,exact_signflip_P=pvalue,interpretation="Exploratory, assumes symmetric donor differences under null; does not address subtype composition or causation"))
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source=s["matrix_path"],count_layer="raw.X",metadata="Use deposited patient/timepoint/material fields; retain original sample IDs. P7 sample aliases map to patient P3, T2 sample aliases to timepoint T1. All primary samples T0. Exact published cell count reproduction unresolved (29 extra cells).",rows=rows,contrasts=contrasts,mainline_validated=False);(B/"results/nmp_qualification_536/PAIRED539.json").write_text(json.dumps(r,indent=2),encoding="utf-8");print(json.dumps(r))
