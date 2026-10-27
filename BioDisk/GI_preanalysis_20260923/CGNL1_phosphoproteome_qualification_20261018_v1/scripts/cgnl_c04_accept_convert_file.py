from pathlib import Path
import sys,subprocess,json
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
name,digest=sys.argv[1:]
r=subprocess.run(['python3',str(P/'scripts/cgnl_c04_accept_convert_only_v1.py'),name,digest],capture_output=True,text=True)
if r.returncode:
 print(r.stdout);print(r.stderr,file=sys.stderr);raise SystemExit(r.returncode)
z=json.loads(r.stdout.strip().splitlines()[-1])
assert z['acquisition_pass']
w=P/'continuous_exploration_20261023_v1/C04_full_fraction_qualification/converted'/Path(name).stem
with (w/'search_pipeline.log').open('w',encoding='utf-8') as log:
 q=subprocess.run(['python3',str(P/'scripts/cgnl_c04_search_full_file.py'),name],stdout=log,stderr=subprocess.STDOUT)
if q.returncode:
 print('Search stage failed; raw/conversion accepted and preserved',file=sys.stderr);raise SystemExit(q.returncode)
z['file_search_complete']=True
print(json.dumps(z))
