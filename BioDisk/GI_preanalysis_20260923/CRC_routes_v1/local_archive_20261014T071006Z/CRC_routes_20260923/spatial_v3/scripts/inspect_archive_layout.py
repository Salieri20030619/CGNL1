from pathlib import Path
import tarfile,sys
sys.stdout.reconfigure(encoding='utf-8')
p=Path(r'C:\Users\31937\CRC_routes_20260923\raw\M-ST-13.tar.gz')
with tarfile.open(p) as t:
 for m in t:print(m.name,m.size,m.offset)
