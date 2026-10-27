from pathlib import Path
import json,datetime
import anndata as ad
import numpy as np
from scipy import sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
state=json.loads((B/'metadata/PROCESSED525_STATE.json').read_text(encoding='utf-8'))
assert state['stage']=='qualified_file_pending_scientific_review'
assert state['actual_sha256']==state['file_manifest']['content_details']['sha256_hash']
spec=dict(purpose='Patient-level assay and lineage qualification; no HGP inference',population='All author-retained nuclei, with author Cell_type labels inspected against marker panels',normalization='CP10k denominator excludes CGNL1',primary_biological_P=None,source_counts='Author-decontaminated expression matrix; verify numeric semantics',endothelial_markers=['PECAM1','CDH5','CLDN5','ESAM','KDR','EMCN','VWF'],sinusoidal_markers=['CLEC4G','STAB2','FCGR2B','LYVE1','CLEC4M','OIT3','DNASE1L3'],no_target_based_selection=True)
with (B/'metadata/INDEPENDENT_QUALIFICATION_SPEC_526.json').open('x',encoding='utf-8') as f:json.dump(spec,f,indent=2)
a=ad.read_h5ad(state['verified_path']);m=sparse.csr_matrix(a.X);m.eliminate_zeros()
assert a.obs_names.is_unique and a.var_names.is_unique and np.isfinite(m.data).all() and (m.data>=0).all()
assert 'CGNL1' in a.var_names
ti=a.var_names.get_loc('CGNL1');target=m[:,ti].toarray().ravel();totals=np.asarray(m.sum(axis=1)).ravel()-target
target_norm=np.divide(10000*target,totals,out=np.zeros(len(target)),where=totals>0)
markers=spec['endothelial_markers']+spec['sinusoidal_markers']+['EPCAM','KRT8','KRT18','KRT19','PTPRC','LST1','RGS5','PDGFRB','ALB']
marker_values={g:m[:,a.var_names.get_loc(g)].toarray().ravel() for g in markers if g in a.var_names}
endo=np.asarray(a.obs['Cell_type'].astype(str)=='Endothelial');patient=np.asarray(a.obs['Patient'].astype(str))
def stats(mask):
    return dict(nuclei=int(mask.sum()),CGNL1_positive_nuclei=int((target[mask]>0).sum()),CGNL1_sum=float(target[mask].sum()),CGNL1_mean_CP10k=float(target_norm[mask].mean()) if mask.any() else None,median_total_excluding_CGNL1=float(np.median(totals[mask])) if mask.any() else None)
by_patient=[]
for p in sorted(set(patient)):
    mask=(patient==p)&endo
    meta=a.obs.loc[patient==p]
    by_patient.append(dict(patient=p,all_nuclei=int((patient==p).sum()),therapy=sorted(set(meta['Therapy'].astype(str))),endothelial=stats(mask),marker_positive_counts={g:int((v[mask]>0).sum()) for g,v in marker_values.items()}))
groups={c:stats(np.asarray(a.obs['Cell_type'].astype(str)==c)) for c in sorted(set(a.obs['Cell_type'].astype(str)))}
subtypes={c:stats(endo&np.asarray(a.obs['Cell_subtype'].astype(str)==c)) for c in sorted(set(a.obs.loc[endo,'Cell_subtype'].astype(str)))}
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),shape=list(m.shape),source_sha256=state['actual_sha256'],all_nonnegative_finite=True,fraction_noninteger_stored_values=float((np.abs(m.data-np.rint(m.data))>1e-6).mean()),independent_patient_labels=len(set(patient)),by_patient=by_patient,by_cell_type=groups,endothelial_subtypes=subtypes,missing_markers=[g for g in markers if g not in marker_values],patients_with_at_least_20_EC=sum(x['endothelial']['nuclei']>=20 for x in by_patient),patients_with_at_least_50_EC=sum(x['endothelial']['nuclei']>=50 for x in by_patient),HGP_field_present=any('hgp' in c.lower() or 'growth' in c.lower() for c in a.obs.columns),inference='Descriptive measurement and lineage qualification only; no HGP association or causal mechanism claim',biological_P=None,mainline_validated=False)
out=B/'results/independent_qualification_526';out.mkdir(exist_ok=False)
(out/'QUALIFICATION526.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='by_patient'}))
print(json.dumps([dict(patient=r['patient'],therapy=r['therapy'],endothelial=r['endothelial']) for r in by_patient]))

