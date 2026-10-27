from pathlib import Path
import json,hashlib,zipfile
root=Path(r'C:\Users\31937\Downloads\单细胞身份_匿名盲审_20261002_v1')
private=Path(r'C:\Users\31937\CGNL1_identity_review_20261002\blind_handoff_private_20261002')
f=root/'00_盲审说明_先读.txt';s=f.read_text(encoding='utf-8');s=s.replace('使用多基因组合及相反证据，不以单个标记或固定“至少两个标记”规则自动决定身份。','使用多基因组合及相反证据，不用单一标记或机械计数门槛代替综合身份判断。');f.write_text(s,encoding='utf-8')
blocked=['CGNL1','CLDN5','ESAM','TJP1','OCLN','GSE245552','GSE285991','GSM','one_marker_removed','two_plus_retained','multiome_prediction','normal_liver_prediction','both_immune','both_endothelial','至少两个标记']
manifest=[]
for f in root.rglob('*'):
 if not f.is_file() or f.name=='MANIFEST.json':continue
 if f.suffix in ['.txt','.tsv','.json']:
  s=f.read_text(encoding='utf-8');assert not any(x in s for x in blocked),f.name
 manifest.append(dict(path=f.relative_to(root).as_posix(),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(root/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=root.with_suffix('.zip')
with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob('*')):
  if f.is_file():z.write(f,arcname=f.relative_to(root))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
qa=json.loads((private/'PACKET_QA.json').read_text(encoding='utf-8'));qa['zip_sha256']=hashlib.sha256(zpath.read_bytes()).hexdigest();qa['specific_previous_filter_rule_hidden']=True;(private/'PACKET_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print('Finalized: 44 clusters, 1320 marker rows, 10 upload files; source labels and previous filter rule hidden.')