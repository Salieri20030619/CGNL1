from pathlib import Path
import hashlib,zipfile,json
r=Path(r'C:\Users\31937\CGNL1_mainline_search_20261012')
p=r/'CGNL1_mainline_search_bundle.zip'
digest=hashlib.sha256(p.read_bytes()).hexdigest()
assert digest=='f8922502be9e98a41fd27af2e60c282be596002d94d8f76347701b330fd938a1'
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 checked=0
 for line in z.read('SHA256SUMS.txt').decode('utf-8').splitlines():
  h,n=line.split(maxsplit=1)
  assert hashlib.sha256(z.read(n)).hexdigest()==h,n
  checked+=1
 for n in ['CGNL1_deep_mainline_review_20261012.html','evidence_register.json','RESEARCH_DECISION.json']:
  assert z.read(n)==(r/n).read_bytes(),n
result={'archive_sha256':digest,'archive_bytes':p.stat().st_size,'manifest_files_verified':checked,'failed':0,'local_report_matches_server':True}
(r/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
