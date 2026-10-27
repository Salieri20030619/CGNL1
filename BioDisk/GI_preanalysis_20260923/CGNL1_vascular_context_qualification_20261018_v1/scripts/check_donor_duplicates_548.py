from pathlib import Path
import anndata as ad,numpy as np,json,hashlib,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");s=json.loads((B/"metadata/ATLAS543_STATE.json").read_text());a=ad.read_h5ad(next(x["path"] for x in s["files"] if x["kind"]=="ec"));m=a.layers["counts"].tocsr();m.sort_indices();groups={}
for donor in ["HD2","HD3"]:
 idx=np.flatnonzero(np.asarray((a.obs.donor==donor)&(a.obs.organ_uni=="liver")));hashes=[]
 for i in idx:
  lo,hi=m.indptr[i:i+2];hashes.append(hashlib.sha256(m.indices[lo:hi].astype("<i8").tobytes()+m.data[lo:hi].astype("<f8").tobytes()).hexdigest())
 groups[donor]=dict(indices=idx.tolist(),row_hashes=hashes,barcodes=a.obs_names[idx].tolist(),sample_ids=sorted(set(a.obs.iloc[idx].sample_id.astype(str))))
x,y=groups["HD2"],groups["HD3"];common=set(x["row_hashes"])&set(y["row_hashes"]);pairs=[]
for h in sorted(common)[:5]:pairs.append(dict(HD2_barcode=x["barcodes"][x["row_hashes"].index(h)],HD3_barcode=y["barcodes"][y["row_hashes"].index(h)],sha256=h))
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),issue="ATLAS-Q01",HD2_n=len(x["indices"]),HD3_n=len(y["indices"]),HD2_sample_ids=x["sample_ids"],HD3_sample_ids=y["sample_ids"],shared_unique_complete_count_vector_hashes=len(common),HD2_unique_hashes=len(set(x["row_hashes"])),HD3_unique_hashes=len(set(y["row_hashes"])),ordered_vectors_identical=x["row_hashes"]==y["row_hashes"],example_pairs=pairs,interpretation="Do not count HD2 and HD3 as independent biological replicates if complete count-vector duplication is confirmed. Source preserved; original attribution requires source clarification.")
(B/"results/normal_atlas_544/DUPLICATE_AUDIT548.json").write_text(json.dumps(r,indent=2),encoding="utf-8");print(json.dumps(r))
