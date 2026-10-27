from pathlib import Path
import json,gzip,csv,hashlib,datetime
O=Path.cwd()/'mechanism_investigation_v2';m=json.loads((O/'metadata/GRCh37_CGNL1_expanded.json').read_text());expected={t['id'] for t in m['Transcript']};files=sorted((O/'raw/GSE237876').glob('GSM765585[12357]_*.txt.gz'));rows=[]
for p in files:
 with gzip.open(p,'rt',encoding='utf-8') as f:
  rd=csv.DictReader(f,delimiter='\t');ids=set();matched=[]
  for r in rd:
   assert r['Name'] not in ids;ids.add(r['Name'])
   if r['Name'] in expected:matched.append({'transcript':r['Name'],'annotated_length':int(r['Length'])})
 rows.append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'total_unique_transcripts':len(ids),'target_rows':matched,'all3targettranscripts_present':len(matched)==3,'target_expression_values_not_used_in_this_check':True})
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':'verifieddownloadedearlyfilesforidentifierqualificationonly','files':rows};(O/'metadata/ROUTE_TARGET_COVERAGE_95.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
