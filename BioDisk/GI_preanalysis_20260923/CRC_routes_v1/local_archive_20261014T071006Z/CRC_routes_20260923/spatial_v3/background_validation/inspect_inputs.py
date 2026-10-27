from pathlib import Path
import gzip,re,sys
sys.stdout.reconfigure(encoding='utf-8')
B=Path(r'C:\Users\31937\CRC_routes_20260923')
for g in ['GSE231559','GSE285991']:
 s=gzip.decompress((B/'metadata'/(g+'.soft.gz')).read_bytes()).decode();print(g)
 for block in s.split('^SAMPLE = ')[1:]:
  title=re.search(r'!Sample_title = (.+)',block).group(1).strip()
  if g=='GSE285991' and not re.fullmatch(r'P\d+_LM',title):continue
  species=re.search(r'!Sample_organism_ch1 = (.+)',block)
  urls=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (.+)',block)
  print(title, species.group(1).strip() if species else '?',urls)
