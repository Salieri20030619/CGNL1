from pathlib import Path
import csv,gzip,json,datetime
import numpy as np
from scipy import io,sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
spec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_508.json').read_text(encoding='utf-8'));q=Path(spec['output'])
with gzip.open(B/'results/pilot_cell_qualification_509/cell_qualification.tsv.gz','rt',encoding='utf-8') as f: rows=list(csv.DictReader(f,delimiter='\t'))
m=sparse.csr_matrix(io.mmread(q/'alevin/quants_mat.mtx'));m.eliminate_zeros()
genes=(q/'alevin/quants_mat_cols.txt').read_text(encoding='utf-8').splitlines()
with (Path(spec['reference']).parent/'transcript_to_gene.tsv').open(encoding='utf-8') as f: lookup={r['gene_id']:r['gene_symbol'] for r in csv.DictReader(f,delimiter='\t')}
symbols=[lookup[g] for g in genes]
total=np.array([float(r['QC_estimated_UMI_excluding_CGNL1']) for r in rows]);det=np.array([int(r['QC_detected_genes_excluding_CGNL1']) for r in rows]);mt=np.array([float(r['mitochondrial_fraction']) for r in rows])
def summary(x): return dict(zip(['min','p25','median','p75','p90','p99','max'],[float(v) for v in np.quantile(x,[0,.25,.5,.75,.9,.99,1])]))
gene_sum=np.asarray(m.sum(axis=0)).ravel();order=np.argsort(gene_sum)[::-1][:20]
markers=json.loads((B/'metadata/LSEC_CONTEXT_REVIEW_SPEC_507.json').read_text(encoding='utf-8'))['diagnostic_markers']
markers+=['CGNL1','PECAM1','CDH5','CLDN5','ESAM','KDR','EMCN','VWF','EPCAM','KRT8','KRT18','KRT19','PTPRC','ALB']
marker_data={}
for g in markers:
    indices=[i for i,s in enumerate(symbols) if s==g]
    if len(indices)!=1: marker_data[g]=dict(mapping_count=len(indices));continue
    values=m[:,indices[0]].toarray().ravel()
    marker_data[g]=dict(mapping_count=1,estimated_UMI_sum=float(values.sum()),positive_cells=int((values>0).sum()))
metrics={name:json.loads((q/name).read_text(encoding='utf-8')) for name in ['aux_info/meta_info.json','aux_info/alevin_meta_info.json','lib_format_counts.json','cmd_info.json'] if (q/name).exists()}
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cells=len(rows),QC_distributions=dict(non_target_estimated_UMI=summary(total),non_target_detected_genes=summary(det),mitochondrial_fraction=summary(mt)),QC_pass_each=dict(UMI500=int((total>=500).sum()),genes200=int((det>=200).sum()),mito20percent=int((mt<=.2).sum()),all=int(((total>=500)&(det>=200)&(mt<=.2)).sum())),top_genes=[dict(gene_id=genes[i],symbol=symbols[i],estimated_UMI_sum=float(gene_sum[i])) for i in order],markers=marker_data,metrics=metrics,primary_QC_unchanged=True,biological_inference=False)
p=B/'results/pilot_numeric_audit_510/QUALITY_DIAGNOSTIC_511.json';p.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

