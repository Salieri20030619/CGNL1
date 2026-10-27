from pathlib import Path
import json,re,hashlib,datetime
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
raw=B/'raw/pilot_source_review_518'
text=(raw/'supplement1.txt').read_text(encoding='utf-8')
table=text.split('Supplementary Table 14.',1)[1].split('Supplementary Table 15.',1)[0]
before,after=table.split('Quality metrics scRNAseq, after manual quality control',1)
def parse(block,fields):
    records=[]
    for line in block.splitlines():
        m=re.match(r'^(PT-\d+)\s+(tumor|adjunct normal)\s+(.+)$',line.strip())
        if not m:continue
        values=m.group(3).split();assert len(values)==len(fields),(line,values)
        row=dict(patient_id=m.group(1),sample_type=m.group(2))
        for key,value in zip(fields,values):row[key]=value if value.endswith('%') else int(value)
        records.append(row)
    assert len(records)==12
    return records
pre=parse(before,['reads','UMIs','events','mean_UMI_per_event','median_UMI_per_event','mean_genes_per_event','median_genes_per_event','saturation'])
post=parse(after,['UMIs','cells','mean_UMI_per_cell','median_UMI_per_cell','mean_genes_per_cell','median_genes_per_cell'])
assert sum(x['cells'] for x in post)==22419
pilot=next(x for x in post if x['patient_id']=='PT-55' and x['sample_type']=='tumor')
assert pilot['cells']==12
tumors=[x for x in post if x['sample_type']=='tumor']
best=max(tumors,key=lambda x:x['median_genes_per_cell'])
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_url='https://media.springernature.com/full/springer-static/esm/art%3A10.1186%2Fs12943-023-01713-1/MediaObjects/12943_2023_1713_MOESM1_ESM.pdf',source_pdf_sha256=hashlib.sha256((raw/'supplement1.pdf').read_bytes()).hexdigest(),page=14,visually_verified=True,before_manual_QC=pre,after_manual_QC=post,post_QC_sum_matches_article=True,PT55_low_recovery_corroborated_by_source=True,tumors_below_50_total_cells=[x['patient_id'] for x in tumors if x['cells']<50],highest_median_gene_recovery_tumor=best['patient_id'],next_pilot_selection_basis='If another pilot is justified, prioritize source median gene recovery and total recoverable cells, not raw download size or CGNL1 expression',tumor_specific_EC_counts_known=False,figureS3_limit='Panel f pools tumor and adjacent-normal cells by patient; it cannot establish tumor-specific EC recovery',decision='Do not use this cohort as sole confirmation of a tumor-EC HGP mechanism. Qualify independent processed data and disease-specific evidence before allocating another raw pilot.',mainline_validated=False)
(B/'metadata/SOURCE_SAMPLE_QUALIFICATION_520.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(pilot=pilot,ranked_tumors=sorted(tumors,key=lambda x:x['median_genes_per_cell'],reverse=True),post_QC_total=22419)))

