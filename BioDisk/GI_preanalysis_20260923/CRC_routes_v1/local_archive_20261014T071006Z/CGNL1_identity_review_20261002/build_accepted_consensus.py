from pathlib import Path
import pandas as pd,json,hashlib,tarfile,datetime
root=Path(r'C:\Users\31937\CGNL1_identity_review_20261002');out=root/'accepted_consensus_20261002_v1';out.mkdir(exist_ok=False)
m=pd.read_csv(root/'external_review_deepseek_20261002_v1/44簇三方解盲比较.tsv',sep='\t',dtype=str,keep_default_na=False)
m['accepted_major_lineage']=m.harmonized_web
m.loc[m.cluster=='C003','accepted_major_lineage']='T_cell'
m['accepted_state']='';m.loc[m.cluster.isin(['C003','C006']),'accepted_state']='proliferating';m['cell_level_caveat']='cluster identity does not certify every barcode'
m.loc[m.cluster=='C038','cell_level_caveat']='epithelial origin; intact-cell and low-RNA QC unresolved'
m[['dataset','cluster','original_dataset','original_cluster','accepted_major_lineage','accepted_state','cell_level_caveat','broad_identity_web','broad_identity_deepseek','confidence_web','confidence_deepseek']].to_csv(out/'accepted_cluster_annotation.tsv',sep='\t',index=False)
receipt={'accepted_at':datetime.datetime.now().astimezone().isoformat(),'user_instruction':'那就按这个取信继续往下','scope':'adopt reviewed major lineages; keep subtype uncertainty and cell-level ambiguity; proceed with fixed target models','first_review_context':'historical context exposure disclosed','second_review_context':'DeepSeek declares no prior exposure; declaration retained','target_retest_not_used_for_annotation':True,'endothelial_clusters':{'GSE245552':'21','GSE285991':'14'},'QC_correction':'GSE245552 mitochondrial metric unavailable with current gene features, not measured zero'}
(out/'ACCEPTANCE.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
files=[(out/'accepted_cluster_annotation.tsv','accepted_cluster_annotation.tsv'),(out/'ACCEPTANCE.json','ACCEPTANCE.json'),(root/'external_review_deepseek_20261002_v1/frozen_original/deepseek_markdown_20261002_56c7cf.md','deepseek_original.md'),(root/'external_review_deepseek_20261002_v1/IMPORT_RECEIPT.json','deepseek_import_receipt.json'),(root/'external_review_context_exposed_20261002_v1/frozen_original/匿名身份审阅_44簇_v1.tsv','first_review_original.tsv'),(root/'external_review_context_exposed_20261002_v1/frozen_original/匿名身份审阅_报告_v1.txt','first_review_original_report.txt'),(root/'external_review_context_exposed_20261002_v1/frozen_original/SHA256_FREEZE.json','first_review_freeze.json')]
manifest=[dict(name=n,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p,n in files];(out/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
with tarfile.open(out/'accepted_consensus.tar.gz','w:gz') as t:
 for p,n in files:t.add(p,arcname=n)
 t.add(out/'MANIFEST.json',arcname='MANIFEST.json')
print(out/'accepted_consensus.tar.gz')