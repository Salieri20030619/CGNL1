from pathlib import Path
import gzip,re,sys
sys.stdout.reconfigure(encoding='utf-8')
B=Path(r'C:\Users\31937\CRC_routes_20260923')
for g in ['GSE231559','GSE285991']:
 s=gzip.decompress((B/'metadata'/(g+'.soft.gz')).read_bytes()).decode()
 print(g)
 for line in s.splitlines():
  if line.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Sample_title','!Sample_characteristics')):print(line[:1000])
