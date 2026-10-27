from pathlib import Path
import gzip,json,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';B=O.parent.parent;out={}
for n in ['GSE314812_sample_design.tsv','GSE314812_paired_design.tsv']:
 p=B/'01_metadata'/n;d=pd.read_csv(p,sep='\t');out[n]=d.to_dict(orient='records');print(n,'shape',d.shape,'columns',list(d));print(d.to_string(index=False)[:10000])
p=B/'01_metadata/GSE314812.soft.gz';b=p.read_bytes();s=(gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b).decode('utf-8');rows=[]
for block in s.split('^SAMPLE = ')[1:]:
 ls=block.splitlines();rows.append({'sample':ls[0],'fields':[l for l in ls if l.startswith(('!Sample_title','!Sample_characteristics_ch1'))]})
out['original_GEO_samples']=rows;(O/'metadata/PAIR_COUNT_AUDIT_104.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print('GEO55_METADATA',json.dumps(rows,ensure_ascii=False)[:12000])
