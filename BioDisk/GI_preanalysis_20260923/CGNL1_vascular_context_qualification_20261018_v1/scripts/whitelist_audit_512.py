from pathlib import Path
import csv,gzip,json,datetime
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
q=Path(json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_508.json').read_text(encoding='utf-8'))['output'])
whitelist=(q/'alevin/whitelist.txt').read_text(encoding='utf-8').splitlines();assert len(set(whitelist))==len(whitelist)
selected=set(whitelist)
with gzip.open(B/'results/pilot_cell_qualification_509/cell_qualification.tsv.gz','rt',encoding='utf-8') as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert selected <= {r['barcode'] for r in rows}
groups={}
for name,rr in [('final_whitelist',[r for r in rows if r['barcode'] in selected]),('matrix_rows_outside_final_whitelist',[r for r in rows if r['barcode'] not in selected])]:
    groups[name]=dict(barcodes=len(rr),QC_pass=sum(r['QC_pass']=='True' for r in rr),endothelial_candidates=sum(r['endothelial_candidate']=='True' for r in rr),ambiguity_filtered_endothelial=sum(r['ambiguity_filtered_endothelial']=='True' for r in rr),CGNL1_estimated_UMI_sum=sum(float(r['CGNL1_estimated_UMI']) for r in rr))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),issue='PILOT-Q02',original_label_correction='1928 represents all emitted matrix barcodes, not final called cells; final whitelist has 1253 barcodes',groups=groups,old_results_preserved=True,biological_negative_claim=False)
(B/'results/pilot_numeric_audit_510/WHITELIST_AUDIT_512.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

