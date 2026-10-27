from pathlib import Path
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');p=R/'scripts'/'finalize_delivery.py';s=p.read_text(encoding='utf-8-sig');block='''# Refresh final input/script checksums after all text/code edits.
manifest=[]
for top in ['01_metadata','02_raw','scripts']:
 for f in sorted((R/top).glob('*')):
  if not f.is_file() or f.name=='integrity_manifest.tsv' or f.name.endswith('.part'):continue
  h=hashlib.sha256()
  with f.open('rb') as stream:
   while True:
    chunk=stream.read(8*1024*1024)
    if not chunk:break
    h.update(chunk)
  manifest.append({'path':str(f.relative_to(R)),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
with (R/'01_metadata'/'integrity_manifest.tsv').open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=['path','bytes','sha256'],delimiter='\\t');writer.writeheader();writer.writerows(manifest)
'''
s=s.replace("pack=D/'GC_validation_evidence_package.zip'",block+"\npack=D/'GC_validation_evidence_package.zip'");p.write_text(s,encoding='utf-8')
