from pathlib import Path
import anndata as ad,numpy as np,json,datetime
from scipy import sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
state=json.loads((B/'metadata/PROCESSED525_STATE.json').read_text(encoding='utf-8'))
a=ad.read_h5ad(state['verified_path']);m=sparse.csr_matrix(a.X,dtype=np.float64);m.eliminate_zeros()
assert np.isfinite(m.data).all() and (m.data>=0).all()
t=m[:,a.var_names.get_loc('CGNL1')].toarray().ravel();denom=np.asarray(m.sum(axis=1)).ravel()-t
norm=np.divide(t*10000,denom,out=np.zeros(len(t)),where=denom>0)
patient=np.asarray(a.obs['Patient'].astype(str));ctype=np.asarray(a.obs['Cell_type'].astype(str));subtype=np.asarray(a.obs['Cell_subtype'].astype(str))
thresholds=[0,1e-8,1e-6,1e-4,.5,1.]
def summarize(mask):
    return dict(nuclei=int(mask.sum()),target_sum=float(t[mask].sum()),mean_CP10k=float(norm[mask].mean()) if mask.any() else None,positive_counts_gt={str(c):int((t[mask]>c).sum()) for c in thresholds})
rows=[]
for p in sorted(set(patient)):
    for s in sorted(set(subtype[ctype=='Endothelial'])):
        mask=(patient==p)&(ctype=='Endothelial')&(subtype==s)
        rows.append(dict(patient=p,subtype=s,**summarize(mask)))
marker_names=['PECAM1','CDH5','CLDN5','ESAM','KDR','EMCN','VWF','PROX1','PDPN','FLT4','LYVE1','CLEC4G','STAB2','FCGR2B','EPCAM','KRT8','PTPRC','RGS5','PDGFRB']
markers={}
for g in marker_names:
    if g not in a.var_names:markers[g]=dict(missing=True);continue
    v=m[:,a.var_names.get_loc(g)].toarray().ravel()
    markers[g]={s:dict(nuclei=int(((ctype=='Endothelial')&(subtype==s)).sum()),positive_gt_1e6=int((v[(ctype=='Endothelial')&(subtype==s)]>1e-6).sum()),estimated_count_sum=float(v[(ctype=='Endothelial')&(subtype==s)].sum())) for s in sorted(set(subtype[ctype=='Endothelial']))}
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),issue='INDEP-Q01',correction='Distinguish nonzero stored values from expression above numerical residuals. Original matrix and preliminary summaries retained; no count rounding.',tolerance_primary=1e-6,tolerance_sensitivity=[1e-8,1e-4],threshold_05_and_1='Additional abundance summaries, not replacement biological detection tests',stored_positive_below_1e6=int(((m.data>0)&(m.data<=1e-6)).sum()),stored_positive_total=int(m.nnz),target_positive_value_quantiles=np.quantile(t[t>0],[0,.01,.1,.5,.9,.99,1]).tolist(),by_cell_type={c:summarize(ctype==c) for c in sorted(set(ctype))},endothelial_by_subtype={s:summarize((ctype=='Endothelial')&(subtype==s)) for s in sorted(set(subtype[ctype=='Endothelial']))},by_patient_subtype=rows,marker_checks=markers,biological_P=None,HGP_test_performed=False,mainline_validated=False)
p=B/'results/independent_qualification_526/NUMERIC_REVIEW527.json';p.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ('marker_checks','by_patient_subtype','by_cell_type')}))
print(json.dumps([r for r in rows if r['nuclei']>0]))

