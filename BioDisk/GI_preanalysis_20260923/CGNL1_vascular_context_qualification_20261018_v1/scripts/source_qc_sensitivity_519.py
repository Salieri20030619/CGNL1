from pathlib import Path
import json,csv,gzip,datetime,importlib.util
import numpy as np
from scipy import sparse,io
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
scope=dict(source='https://pmc.ncbi.nlm.nih.gov/articles/PMC9872436/',change='Separate source-threshold sensitivity on existing Alevin matrix; primary analysis and thresholds retained',source_thresholds=dict(min_detected_genes=150,max_mitochondrial_fraction=.7,no_reported_total_UMI_threshold=True,gene_min_cells=3),not_exact_reproduction='Original used CellRanger, graph clustering and iterative manual removal; this only tests reported initial QC thresholds with our frozen identity rules',biological_P=None,post_primary_inspection=True,not_a_replacement_primary_gate=True)
with (B/'metadata/SOURCE_QC_SENSITIVITY_SPEC_519.json').open('x',encoding='utf-8') as f:json.dump(scope,f,indent=2)
q=Path(json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_508.json').read_text(encoding='utf-8'))['output'])
cells=(q/'alevin/quants_mat_rows.txt').read_text(encoding='utf-8').splitlines();genes=(q/'alevin/quants_mat_cols.txt').read_text(encoding='utf-8').splitlines()
mat=sparse.csr_matrix(io.mmread(q/'alevin/quants_mat.mtx'));mat.eliminate_zeros()
spec=json.loads((B/'metadata/PILOT_CELL_QUALIFICATION_AMENDMENT_458.json').read_text(encoding='utf-8'))
qspec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_508.json').read_text(encoding='utf-8'))
with (Path(qspec['reference']).parent/'transcript_to_gene.tsv').open(encoding='utf-8') as f: lookup={r['gene_id']:r['gene_symbol'] for r in csv.DictReader(f,delimiter='\t')}
symbols=[lookup[g] for g in genes]
module_spec=importlib.util.spec_from_file_location('frozen_classifier',B/'scripts/qualify_pilot_cells_509.py');module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
spec['cell_QC_primary']=dict(min_detected_genes=150,max_mitochondrial_UMI_fraction=.7,min_total_UMI=0)
res=module.classify(mat,symbols,spec)
wl=set((q/'alevin/whitelist.txt').read_text(encoding='utf-8').splitlines());called=np.array([c in wl for c in cells])
gene_support=np.asarray((mat[called,:]>0).sum(axis=0)).ravel()>=3
gene_support[np.array(symbols)=='CGNL1']=False
det_filtered=np.asarray((mat[:,gene_support]>0).sum(axis=1)).ravel()
qc_filtered=res['qc']&(det_filtered>=150)
ti=symbols.index('CGNL1');target=mat[:,ti].toarray().ravel()
groups={}
for name,mask in [('final_whitelist',called),('other_matrix_barcodes',~called)]:
    groups[name]=dict(barcodes=int(mask.sum()),initial_source_threshold_QC=int((mask&res['qc']).sum()),source_QC_after_gene_min3=int((mask&qc_filtered).sum()),EC_candidates_source_threshold=int((mask&res['candidate']).sum()),EC_candidates_after_gene_min3=int((mask&qc_filtered&res['candidate']).sum()),ambiguity_filtered_EC_source_threshold=int((mask&res['qualified']).sum()),CGNL1_estimated_UMI_sum=float(target[mask].sum()))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=scope,groups=groups,primary_unchanged=True,mainline_validated=False)
(B/'results/pilot_numeric_audit_510/SOURCE_QC_SENSITIVITY_519.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

