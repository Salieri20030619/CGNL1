from pathlib import Path
import json,subprocess
O=Path.cwd()/'mechanism_investigation_v2';logs=sorted((O.parent/'logs').glob('cgnl_mech177_server_failed_file_recovery*.log'));print('RECOVERY_LOG',logs[-1].read_text(encoding='utf-8')[-3500:] if logs else 'missing');print('RECOVERY_EXIT_FILES',[p.name for p in (O.parent/'logs').glob('cgnl_mech177*.exit.json')]);print('PROXY_LISTENER',[l for l in subprocess.run(['ss','-ltnp'],capture_output=True,text=True).stdout.splitlines() if ':7897' in l]);print('RECOVERY_PROC',Path('/proc/869788').exists())
