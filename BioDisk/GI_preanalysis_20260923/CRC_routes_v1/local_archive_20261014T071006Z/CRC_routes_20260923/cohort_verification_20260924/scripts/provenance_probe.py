from pathlib import Path
import gzip,json
B=Path(r'C:\Users\31937\CRC_routes_20260923'); A=B/'cohort_verification_20260924';out=[]
for g in ['GSE225857','GSE231559','GSE285991','GSE245552','GSE178318']:
 p=B/'metadata'/(g+'.soft.gz')
 if not p.exists():p=A/'sc_metadata'/(g+'.soft.gz')
 rec={'dataset':g,'metadata':{},'GSM':[]}
 with gzip.open(p,'rt',errors='replace') as f:
  for line in f:
   if line.startswith(('!Series_title','!Series_contact_institute','!Series_contact_name','!Series_pubmed_id','!Series_submission_date','!Series_sample_id')):
    k,v=line.strip().split(' = ',1);rec['metadata'].setdefault(k,[]).append(v)
   if line.startswith('^SAMPLE = '):rec['GSM'].append(line.strip().split(' = ')[1])
 out.append(rec)
sets=[set(x['GSM']) for x in out]; overlaps=[]
for i in range(len(out)):
 for j in range(i):
  ov=sets[i]&sets[j]
  if ov:overlaps.append({'a':out[i]['dataset'],'b':out[j]['dataset'],'GSM':sorted(ov)})
res={'studies':out,'cross_study_GSM_overlap':overlaps,'boundary':'Distinct public library accessions and donor labels do not provide global person-level identity linkage.'};(A/'SC_PROVENANCE.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
