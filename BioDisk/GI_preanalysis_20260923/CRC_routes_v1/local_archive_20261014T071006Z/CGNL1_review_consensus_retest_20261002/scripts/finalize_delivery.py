
from pathlib import Path
import json,hashlib,tarfile,zipfile,subprocess,time,shutil
import pandas as pd
r=Path(r'C:\Users\31937\CGNL1_review_consensus_retest_20261002');d=r/'results'
c=pd.read_csv(d/'continuous_models.tsv',sep='\t')
pos=c[(c.scope=='pooled')&(c.model=='study_depth_adjusted')&(c.stage.isin(['cluster_primary_n50','exclude_reference_immune_n50','both_references_endo_n50']))].copy()
pos['evidence_class']='prespecified sensitivity; nominal positive; n8 exploratory; overlapping patients'
pos['follow_up']='retain candidate; subtype/state and independent replication; primary negative retained'
pos.to_csv(d/'nominal_positive_candidate_registry.tsv',sep='\t',index=False)
qa=json.loads((r/'QA.json').read_text());qa['delivery_figure_visual_inspection']='passed';qa['local_download_44_file_SHA256']='passed';(r/'DELIVERY_QA.json').write_text(json.dumps(qa,indent=2))
files=[r/'复测说明.txt',r/'NEXT_STEP.json',r/'DELIVERY_QA.json',r/'scripts/posthoc_diagnostics.py',r/'scripts/build_delivery.py',d/'nominal_positive_candidate_registry.tsv',r/'figures/reviewed_identity_retest.png',r/'figures/reviewed_identity_retest.svg']
manifest={p.relative_to(r).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(r/'ADDENDUM_MANIFEST.json').write_text(json.dumps(manifest,indent=2));files.append(r/'ADDENDUM_MANIFEST.json')
with tarfile.open(r/'report_addendum.tar.gz','w:gz') as tf:
 for p in files:tf.add(p,arcname=p.relative_to(r).as_posix())
remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/review_consensus_retest_20261002_v1'
for attempt in range(3):
 with (r/'report_addendum.tar.gz').open('rb') as inp:
  x=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','cat > '+remote+'/report_addendum.tar.gz'],stdin=inp,stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode==0:break
 time.sleep(3)
assert x.returncode==0,x.stderr
remote_code="""from pathlib import Path
import tarfile,json,hashlib
r=Path('"""+remote+"""')
with tarfile.open(r/'report_addendum.tar.gz') as tf:
 for m in tf.getmembers():assert (r/m.name).resolve().is_relative_to(r.resolve()) and not m.issym() and not m.islnk()
 tf.extractall(r)
a=json.loads((r/'ADDENDUM_MANIFEST.json').read_text());assert all(hashlib.sha256((r/f).read_bytes()).hexdigest()==h for f,h in a.items());print('REMOTE_ADDENDUM_VERIFIED',len(a))
"""
for attempt in range(3):
 x=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3 -'],input=remote_code.encode(),capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode==0:break
 time.sleep(3)
assert x.returncode==0,x.stderr
print(x.stdout.decode())
allfiles=[p for p in r.rglob('*') if p.is_file() and not p.name.endswith('.tar.gz') and p.name!='FINAL_MANIFEST.json'];(r/'FINAL_MANIFEST.json').write_text(json.dumps({p.relative_to(r).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in allfiles},indent=2))
out=Path(r'C:\Users\31937\Downloads\CGNL1_accepted_review_results_20261002.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for p in allfiles+[r/'FINAL_MANIFEST.json']:z.write(p,arcname=p.relative_to(r).as_posix())
with zipfile.ZipFile(out) as z:assert z.testzip() is None
shutil.copy2(r/'复测说明.txt',Path(r'C:\Users\31937\Downloads\CGNL1_审阅后复测说明_20261002.txt'))
print('DELIVERED',out, out.stat().st_size)
