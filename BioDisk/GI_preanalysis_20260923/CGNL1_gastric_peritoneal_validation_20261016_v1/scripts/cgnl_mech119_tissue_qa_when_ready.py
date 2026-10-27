from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,datetime,tarfile,gzip,hashlib,csv,io
import numpy as np,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';done=O/'results/GSE308231_SOURCE_COMPLETE.json'
def live():
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if any(a.endswith(b'/cgnl_mech47b_tissue_when_ready.py') for a in (p/'cmdline').read_bytes().split(b'\0')):return True
  except (FileNotFoundError,PermissionError):pass
 return False
while not done.exists():
 if not live():raise RuntimeError('Source controller terminal without completion marker; inspect recorded exit')
 time.sleep(900)
receipt=json.loads((O/'metadata/GSE308231_DOWNLOAD.json').read_text(encoding='utf-8'));assert receipt['success'];p=O/'raw/GSE308231_RAW.tar';h=hashlib.sha256()
with p.open('rb') as f:
 while b:=f.read(4194304):h.update(b)
assert h.hexdigest()==receipt['sha256'] and p.stat().st_size==receipt['bytes']
features=[]
with tarfile.open(p) as tf:
 for m in tf.getmembers():
  if not m.name.endswith('_features.tsv.gz'):continue
  b=gzip.decompress(tf.extractfile(m).read());rr=list(csv.reader(io.StringIO(b.decode('utf-8')),delimiter='\t'));ids=[r[0] for r in rr];target=[r for r in rr if r[1]=='CGNL1'];assert len(ids)==len(set(ids));assert len(target)==1 and target[0][0]=='ENSG00000128849'
  features.append({'member':m.name,'rows':len(rr),'human_ENSG':sum(x.startswith('ENSG') for x in ids),'pig_ENSSSCG':sum(x.startswith('ENSSSCG') for x in ids),'target':target,'decompressed_sha256':hashlib.sha256(b).hexdigest()})
assert len(features)==6 and all(x['human_ENSG']==x['rows'] and x['pig_ENSSSCG']==0 for x in features)
a=pd.read_csv(O/'results/GSE308231_source_profiles.tsv',sep='\t');qc=pd.read_csv(O/'results/GSE308231_QC.tsv',sep='\t');assert len(qc)==12 and a['sample'].nunique()==6;assert not a[['sample','scheme','group']].duplicated().any();checks=[]
for row in qc.itertuples():
 q=a[(a['sample']==row.sample)&(a.scheme==row.scheme)];assert len(q)==8
 assert q.n_cells.sum()==row.retained_cells and q.target_UMI.sum()==row.target_UMI
 assert (q.target_detected_cells<=q.n_cells).all() and (q.target_detected_cells<=q.target_UMI).all()
 for r in q.itertuples():
  if r.n_cells==0:assert r.target_UMI==0 and r.total_UMI==0 and pd.isna(r.target_detection)
  else:assert abs(r.target_detection-r.target_detected_cells/r.n_cells)<1e-12
  if r.total_UMI>0:assert abs(r.CGNL1_ln1pCP10k-np.log1p(r.target_UMI/r.total_UMI*1e4))<1e-12
 checks.append({'sample':row.sample,'scheme':row.scheme,'retained_cells':row.retained_cells,'target_UMI':row.target_UMI,'source_partition_accounting_pass':True})
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'archive_sha256':h.hexdigest(),'feature_checks':features,'source_accounting_checks':checks,'scope':'wholearchiveidentity; allsix feature species/target annotation; sample x source accounting and reported transforms','not_verified':['independent biological donor pairing','marker-defined identity by genotype or protein','full replication of original doublet removal'],'metadata_issue_GC_Q07':'human-reference qualification does not resolve inconsistent PEDV processing description'}
(O/'results/GSE308231_SOURCE_QA_119.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
