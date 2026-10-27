from pathlib import Path
import json,sys
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');R=Path(json.loads((B/'LATEST_RUN.json').read_text())['local_run'])
for name in ['v01_strict_masks.py','task_entry.py','data_prepare.py','launch_v01.py']:
 s=(B/name).read_text(encoding='utf-8-sig');compile(s,name,'exec');(R/'scripts'/name).write_text(s,encoding='utf-8')
print('4 scripts syntax-checked')
