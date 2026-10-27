import json,gzip,urllib.request
from pathlib import Path
r=Path(r'C:/Users/31937/CGNL1_data_prepare_20260929');m=Path(r'C:/Users/31937/CRC_routes_20260923/expanded_screen_20260925/metadata')
jobs=[]
for l in gzip.decompress((m/'GSE131418.soft.gz').read_bytes()).decode().splitlines():
 if l.startswith('!Series_supplementary_file'):
  u=l.split(' = ',1)[1].replace('ftp://','https://')
  if not u.endswith('_RAW.tar'):jobs.append({'dataset':'GSE131418','name':u.rsplit('/',1)[-1],'url':u,'role':'bulk_processed_expression_or_clinical_or_probe_annotation'})
x=json.loads((r/'osf_root.json').read_text())['data'][0]
jobs.append({'dataset':'oncoPredict_GDSC_CTRP','name':'DataFiles.zip','url':x['links']['download'],'sha256':x['attributes']['extra']['hashes']['sha256'],'expected_bytes':x['attributes']['size'],'role':'published_drug_training_resources'})
for name in ['README_eQTL_v8.txt','GTEx_Analysis_v8_eQTL.tar','GTEx_Analysis_v8_eQTL_independent.tar']:
 jobs.append({'dataset':'GTEx_v8','name':name,'url':'https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/'+name,'role':'MR_exposure_qualification_not_complete_colocalization'})
for name in ['finngen_R12_C3_COLORECTAL_EXALLC.gz','finngen_R12_C3_COLORECTAL_EXALLC.gz.tbi']:
 jobs.append({'dataset':'FinnGen_R12_CRC','name':name,'url':'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/'+name,'role':'general_CRC_susceptibility_not_liver_metastasis'})
jobs.append({'dataset':'FinnGen_R12_CRC','name':'finngen_R12_pheno_n.tsv','url':'https://storage.googleapis.com/finngen-public-data-r12/finngen_R12_pheno_n.tsv','role':'phenotype_definitions'})
(r/'wave02_extended.json').write_text(json.dumps(jobs,indent=2));print('WAVE2_FILES',len(jobs))
