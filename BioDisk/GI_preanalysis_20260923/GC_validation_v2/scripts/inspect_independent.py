import pathlib,tarfile,json,gzip,pandas as pd,xml.etree.ElementTree as ET
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2';p=R/'02_raw'/'GSE228598_RAW.tar'
with tarfile.open(p) as t:
 members=t.getmembers();print('MEMBERS',len(members));print('\n'.join(f'{x.name}\t{x.size}' for x in members[:15]));(R/'01_metadata'/'GSE228598_archive_members.json').write_text(json.dumps([{'name':x.name,'size':x.size} for x in members],indent=2))
 for mem in members[:3]:
  if mem.isfile():
   f=t.extractfile(mem)
   if mem.name.endswith('.gz'):f=gzip.GzipFile(fileobj=f)
   print(mem.name,[f.readline().decode().strip()[:500] for _ in range(3)])
p=R/'01_metadata'/'GSE162213.soft.gz'
print((R/'01_metadata'/'GSE162213.files.html').read_text())
for l in gzip.open(p,'rt'):
 if l.startswith(('!Series_title','!Series_overall_design','!Sample_title','!Sample_supplementary_file')):print(l.strip())
