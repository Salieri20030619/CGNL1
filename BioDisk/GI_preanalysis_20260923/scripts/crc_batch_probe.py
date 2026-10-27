import pathlib
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
s=None
for l in (r/'01_metadata'/'GSE145356.soft.txt').read_text().splitlines():
 if l.startswith('^SAMPLE'):print(l)
 if l.startswith(('!Sample_title','!Sample_instrument','!Sample_library','!Sample_source')):print(l)
