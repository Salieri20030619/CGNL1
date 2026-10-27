from pathlib import Path
import json,hashlib,zipfile
A=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/cohort_verification_20260924')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
s=json.loads((A/'AUDIT_COMPLETE.json').read_text());z=A/'cohort_verification_evidence_20260924.zip';assert sha(z)==s['evidence_zip_sha256']
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for x in f.namelist():assert (A/x).resolve().is_relative_to(A.resolve())
 f.extractall(A)
for r in json.loads((A/'AUDIT_DELIVERY_MANIFEST.json').read_text()):
 p=A/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['path']
for name,expected in [('GSE159216_series_matrix.txt.gz','79a73b2ac8ceeaa2d5a1d0bf7a631273dc99357fb278d32f77b57717e578c93c'),('GSE73255_series_matrix.txt.gz','2161c93539c1bf22ce5de5a6124fc71c46a2cf6e58915640f843fa69c60a742d')]:assert sha(A/'bulk'/name)==expected
result={'exit_code':0,'evidence_archive_sha256_verified':True,'evidence_files_verified':True,'bulk_matrices_sha256_verified':True,'audit_complete':True}
(A/'SERVER_ARCHIVE_VERIFIED.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
