from common import *
import shutil
D=Q/'V04_INPUT_QC';D.mkdir(exist_ok=True);old=pathlib.Path(C['previous_queue'])/'V04_INPUT_QC';r=json.loads((old/'RESULT.json').read_text());assert r['status']=='whitelist_matrices_ready'
for p in old.iterdir():
 if p.is_file() and p.name not in ['RESULT.json','EXIT.json','COMPLETE.json','STARTED.json','run.log']:
  dest=D/p.name
  if not dest.exists():dest.symlink_to(p)
dump(D/'RESULT.json',{'status':'whitelist_matrices_reused','original':str(old),'samples':2,'clinical_and_morphology_review':'pending'})
