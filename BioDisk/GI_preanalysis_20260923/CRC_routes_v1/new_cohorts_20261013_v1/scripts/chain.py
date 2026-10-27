from pathlib import Path
import subprocess,sys
R=Path(__file__).resolve().parent.parent
for name in ['04_extract_primary.py','05_analyze.py','06_qa.py','07_report.py']:
 p=subprocess.run([sys.executable,str(R/'scripts/run.py'),name]);assert p.returncode==0,name