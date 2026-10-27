from pathlib import Path
import pandas as pd,json
O=Path.cwd()/'mechanism_investigation_v2';a=pd.read_csv(O/'results/PMGC2026_scRNA_source_profiles.tsv',sep='\t');print(a[a.scheme=='primary'].to_string(index=False));print('Sensitivity mesothelial',a[(a.scheme=='sensitivity')&a.group.eq('strict_mesothelial')].to_string(index=False))
q=json.loads((O/'results/TISSUE_RAW_COUNT_QA_151b.json').read_text());n=sum(x['triplets'] for x in q['samples']);ledger=json.loads((O/'results/MECHANISM_EVIDENCE_LEDGER_v9.json').read_text(encoding='utf-8'))
for r in ledger:
 if r.get('id')=='MEC-15':r['quality']='all18gzCRC;6featurefileshumanENSG/targetverified;sourceaccountingQA119passed; independent line-by-line '+str(n)+' MatrixMarket triplets checked in audit151b';r['record_correction']='previousv9 shorthand173million was arithmetic transcription error; actual summed auditlog count '+str(n)+'; expression results unchanged'
(O/'results/MECHANISM_EVIDENCE_LEDGER_v10.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2),encoding='utf-8');p=O/'STATUS.json';s=json.loads(p.read_text(encoding='utf-8'));s['mechanism_ledger']='results/MECHANISM_EVIDENCE_LEDGER_v10.json';p.write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8');print('Correct audit triplets',n)
