from pathlib import Path
import anndata as ad,h5py,numpy as np,json,math,datetime,itertools
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");s=json.loads((B/"metadata/NMP534_STATE.json").read_text());old=json.loads((B/"results/nmp_qualification_536/PAIRED539.json").read_text());a=ad.read_h5ad(s["matrix_path"],backed="r");genes=list(a.raw.var_names);markers=["CGNL1","PECAM1","CDH5","CLDN5","ESAM","PROX1","PDPN","LYVE1","FLT4","CLEC4G","STAB2","FCGR2B","ACKR1","SELP","RGCC","KDR","EPCAM","PTPRC"]
assert len(genes)==len(set(genes));gi={g:genes.index(g) for g in markers};idx=np.flatnonzero(np.asarray((a.obs.cell_type=="Endothelial cells")&(a.obs.timepoint=="T0")));groups={};cells=[]
with h5py.File(s["matrix_path"],"r") as h:
 raw=h["raw/X"];ptr=raw["indptr"][:]
 for i in idx:
  lo,hi=int(ptr[i]),int(ptr[i+1]);indices=raw["indices"][lo:hi];values=raw["data"][lo:hi].astype(float);assert len(indices)==len(set(indices.tolist()));d=dict(zip(indices.tolist(),values.tolist()));v={g:d.get(j,0.) for g,j in gi.items()};obs=a.obs.iloc[i];key=(str(obs.patient),str(obs.material));g=groups.setdefault(key,dict(n=0,target=0.,denominator=0.,positive=0));g["n"]+=1;g["target"]+=v["CGNL1"];g["denominator"]+=math.fsum(values.tolist())-v["CGNL1"];g["positive"]+=int(v["CGNL1"]>0);cells.append(dict(patient=key[0],tissue=key[1],barcode=str(a.obs_names[i]),markers=v))
errors=[]
for row in old["rows"]:
 g=groups[row["patient"],row["tissue"]];effect=math.log2(1+1e6*g["target"]/g["denominator"])
 for k,v in [("n_cells",g["n"]),("CGNL1_count",g["target"]),("non_target_UMI",g["denominator"]),("positive_cells",g["positive"]),("log2CPM1",effect)]:
  if not math.isclose(row[k],v,rel_tol=1e-12,abs_tol=1e-12):errors.append(dict(patient=row["patient"],field=k,old=row[k],recomputed=v))
for contrast in old["contrasts"]:
 d=[x["delta"] for x in contrast["pairs"]];obs=abs(math.fsum(d));tail=sum(abs(math.fsum(x if (mask>>i)&1 else -x for i,x in enumerate(d)))>=obs-1e-12 for mask in range(2**len(d)));p=tail/2**len(d)
 if p!=contrast["exact_signflip_P"]:errors.append(dict(test="signflip",cutoff=contrast["cutoff"]))
summary=[]
for (patient,tissue),g in groups.items():
 cc=[x for x in cells if x["patient"]==patient and x["tissue"]==tissue];positive=[x for x in cc if x["markers"]["CGNL1"]>0]
 summary.append(dict(patient=patient,tissue=tissue,n=len(cc),target_positive=len(positive),marker_positive_counts={m:dict(all=sum(x["markers"][m]>0 for x in cc),target_positive=sum(x["markers"][m]>0 for x in positive)) for m in markers if m!="CGNL1"}))
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),method="Independent direct h5py CSR row reads and Python math.fsum; original computation used anndata/scipy sparse aggregation. Exact test recalculated by bit enumeration.",errors=errors,passed=not errors,marker_summary=summary,interpretation="Descriptive co-detection only; neither absence under dropout nor presence of one marker defines cell identity. No new inferential tests.")
out=B/"results/nmp_qualification_536";(out/"AUDIT540.json").write_text(json.dumps(r,indent=2),encoding="utf-8");(out/"ENDOTHELIAL_MARKERS540.json").write_text(json.dumps(cells,indent=2),encoding="utf-8");print(json.dumps(r));a.file.close()
