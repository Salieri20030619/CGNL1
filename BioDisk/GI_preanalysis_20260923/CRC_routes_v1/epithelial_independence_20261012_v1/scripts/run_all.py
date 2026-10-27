from pathlib import Path
import subprocess,sys,json,datetime
R=Path(__file__).resolve().parent.parent
jobs=[('extract','01_extract.py'),('AB','02_analyze_AB.py'),('B_subtype','02b_verify_subtype.py'),('C','03_analyze_C.py'),('C_adapted','03b_spatial_platform_adapted.py'),('A_mixture','02c_lineage_mixture.py'),('C_blocks','03c_block_diagnostic.py'),('QA','04_numeric_qa.py'),('report','05_report.py')]
status={}
for name,script in jobs:
 with (R/'logs'/f'{name}.log').open('w') as log:rc=subprocess.call([sys.executable,str(R/'scripts'/script)],cwd=R,stdout=log,stderr=subprocess.STDOUT)
 (R/'logs'/f'{name}.exit').write_text(str(rc)+'\n');status[name]=rc
 print(name,rc,flush=True)
 if rc:break
(R/'REPRODUCTION_STATUS.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_codes':status,'all_success':len(status)==len(jobs) and all(x==0 for x in status.values())},indent=2))
sys.exit(0 if len(status)==len(jobs) and all(x==0 for x in status.values()) else 1)
