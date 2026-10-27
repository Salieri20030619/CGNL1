from pathlib import Path
import csv, gzip, hashlib, json, datetime
import numpy as np
from scipy import io, sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
assert json.loads((B/'metadata/CELL509.exit.json').read_text(encoding='utf-8'))['exit_code']==0
spec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_508.json').read_text(encoding='utf-8'))
q=Path(spec['output'])
results=B/'results/pilot_cell_qualification_509'
reported=json.loads((results/'CELL_QUALIFICATION_COMPLETE_509.json').read_text(encoding='utf-8'))
cells=(q/'alevin/quants_mat_rows.txt').read_text(encoding='utf-8').splitlines()
genes=(q/'alevin/quants_mat_cols.txt').read_text(encoding='utf-8').splitlines()
mat=sparse.csr_matrix(io.mmread(q/'alevin/quants_mat.mtx'))
assert mat.shape==(len(cells),len(genes))
assert len(set(cells))==len(cells) and len(set(genes))==len(genes)
assert np.isfinite(mat.data).all() and (mat.data>=0).all()
with (Path(spec['reference']).parent/'transcript_to_gene.tsv').open(encoding='utf-8') as f:
    mapping={r['gene_id']:r['gene_symbol'] for r in csv.DictReader(f,delimiter='\t')}
symbols=[mapping[g] for g in genes]
assert symbols.count('CGNL1')==1
ti=symbols.index('CGNL1')
target=mat[:,ti].toarray().ravel()
non_target=np.array([i for i,s in enumerate(symbols) if s!='CGNL1'])
expected_total=np.asarray(mat[:,non_target].sum(axis=1)).ravel()
expected_detected=np.asarray((mat[:,non_target]>0).sum(axis=1)).ravel()
with gzip.open(results/'cell_qualification.tsv.gz','rt',encoding='utf-8') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
assert [r['barcode'] for r in rows]==cells
assert np.allclose([float(r['CGNL1_estimated_UMI']) for r in rows],target,rtol=1e-12,atol=1e-9)
assert np.allclose([float(r['QC_estimated_UMI_excluding_CGNL1']) for r in rows],expected_total,rtol=1e-12,atol=1e-9)
assert np.array_equal([int(r['QC_detected_genes_excluding_CGNL1']) for r in rows],expected_detected)
flags={'all_called':np.ones(len(cells),dtype=bool)}
for group,column in [('QC_pass','QC_pass'),('endothelial_candidate','endothelial_candidate'),('ambiguity_filtered_endothelial','ambiguity_filtered_endothelial')]:
    assert all(r[column] in ('True','False') for r in rows)
    flags[group]=np.array([r[column]=='True' for r in rows])
for group,mask in flags.items():
    actual=reported['groups'][group]
    assert int(mask.sum())==actual['cells']
    assert int((target[mask]>0).sum())==actual['CGNL1_positive_cells']
    assert np.isclose(target[mask].sum(),actual['CGNL1_estimated_UMI_sum'],rtol=1e-12,atol=1e-9)
    fraction=float((target[mask]>0).mean()) if mask.any() else None
    assert fraction==actual['CGNL1_detection_fraction']
markers=['CGNL1','PECAM1','CDH5','CLDN5','ESAM','KDR','EMCN','VWF','EPCAM','KRT8','KRT18','KRT19','RGS5','PDGFRB','CSPG4','ACTA2','ALB','APOA1','TTR','ASGR1','PTPRC','LST1','TYROBP','CD3D']
marker_summaries={}
for marker in markers:
    assert symbols.count(marker)==1
    values=mat[:,symbols.index(marker)].toarray().ravel()
    marker_summaries[marker]={group:dict(cells=int(mask.sum()),positive_cells=int((values[mask]>0).sum()),estimated_UMI_sum=float(values[mask].sum())) for group,mask in flags.items()}
metrics={}
for name in ['aux_info/meta_info.json','aux_info/alevin_meta_info.json','cmd_info.json']:
    path=q/name
    if path.exists(): metrics[name]=json.loads(path.read_text(encoding='utf-8'))
digest=hashlib.sha256()
with (q/'alevin/quants_mat.mtx').open('rb') as f:
    for block in iter(lambda:f.read(8*1024*1024),b''): digest.update(block)
assert digest.hexdigest()==reported['source_matrix_sha256']
audit=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),numeric_checks_passed=True,independent_donors=1,source_matrix_sha256=digest.hexdigest(),marker_summaries=marker_summaries,quantification_metrics=metrics,manual_identity_review_pending=True,mainline_validated=False)
destination=B/'results/pilot_numeric_audit_510'
destination.mkdir(exist_ok=False)
(destination/'AUDIT510.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(numeric_checks_passed=True,cells=len(cells),genes=len(genes),output=str(destination))))

