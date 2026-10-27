import pathlib,numpy as np,pandas as pd
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');a=np.load(B/'results/route1/aggregates.npz');m=pd.read_csv(B/'results/route1/group_metadata.tsv',sep='\t');sy=a['genes'];L=a['log2CPM'];rows=[]
for gene in ['CGNL1','CREB3L1','PCSK5','SLC9B2','TMEM131L','BTNL8']:
 i=np.flatnonzero(sy==gene)[0]
 for ct in ['Tumor_author','Tumor_strict_QC','Fibroblast','Endothelial']:
  z=m[(m.celltype==ct)&(m.n_cells>=20)].copy();z['value']=L[i,z.index];v=z.groupby('patient').value.mean();rows.append({'gene':gene,'celltype':ct,'n':len(v),'patient_median_log2CPM':v.median()})
r=pd.DataFrame(rows);(B/'deep_v2/results/mechanism_candidates').mkdir(exist_ok=True);r.to_csv(B/'deep_v2/results/mechanism_candidates/source_audit.tsv',sep='\t',index=False);print(r.to_string(index=False))
