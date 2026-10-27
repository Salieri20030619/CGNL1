from pathlib import Path
import json
O=Path.cwd()/'mechanism_investigation_v2'
print('CURRENT_STATUS', (O/'STATUS.json').read_text(encoding='utf-8')[:10000])
print('LEDGER',(O/'results/MECHANISM_EVIDENCE_LEDGER_v8.json').read_text(encoding='utf-8'))
print('RESULT_FILES',json.dumps([{'root':str(p),'files':[x.name for x in p.iterdir() if x.is_file()]} for p in [O.parent/'results',O/'results']],ensure_ascii=False))
