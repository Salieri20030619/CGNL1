from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,zipfile
A=Path(__file__).resolve().parent.parent
old=pd.read_csv(A/'existing_endothelial_coverage.tsv',sep='\t'); old=old.rename(columns={'frozen_marker_endo':'endo_cells','CGNL1_positive_frozen':'CGNL1_positive'})
parts=[old[['dataset','donor','endo_cells','CGNL1_positive','matched_pairs','eligible_pseudobulk50','eligible_matched_state']]]
for g in ['GSE245552','GSE178318']:
 z=pd.read_csv(A/(g+'_library_coverage.tsv'),sep='\t');z=z[z.eligible_pseudobulk50].sort_values(['donor','endo_cells'],ascending=[True,False]).drop_duplicates('donor')
 parts.append(z[parts[0].columns])
q=pd.concat(parts,ignore_index=True);assert len(q)==14 and q.donor.is_unique and q.eligible_matched_state.sum()==6
q.to_csv(A/'qualified_endothelial_patients.tsv',sep='\t',index=False)
rows=[]
for g,n in [('GSE225857',5),('GSE231559',8),('GSE285991',10),('GSE245552',14),('GSE178318',6)]:
 z=q[q.dataset==g]; rows.append({'dataset':g,'LM_patient_labels':n,'eligible_endo_pseudobulk_patients':len(z),'eligible_CGNL1_state_patients':int(z.eligible_matched_state.sum()),'role':'discovery' if g=='GSE225857' else 'external'})
cs=pd.DataFrame(rows);cs.to_csv(A/'cohort_summary.tsv',sep='\t',index=False);assert cs.LM_patient_labels.sum()==43
bulk=[]
for g in ['GSE159216','GSE73255']:
 m=pd.read_csv(A/'bulk'/(g+'_sample_metadata_CGNL1.tsv'),sep='\t',dtype=str,keep_default_na=False)
 if g=='GSE159216':
  z=m;idcol='patient number';tissue='tissue';expected=(283,171)
 else:
  z=m[m['tissue type']=='frozen liver metastasis sample'];idcol='tissue id';tissue='tissue type';expected=(185,167)
 assert (len(z),z[idcol].nunique())==expected and not z[idcol].isin(['','NA']).any()
 cg=[c for c in m.columns if c.startswith('CGNL1__')];finite={c:int(np.isfinite(pd.to_numeric(z[c],errors='coerce')).sum()) for c in cg};assert all(v==len(z) for v in finite.values())
 clinical={}
 for c in [c for c in z.columns if 'survival' in c]:
  valid=pd.to_numeric(z[c],errors='coerce');good=valid.notna();clinical[c]={'patients':int(z.loc[good,idcol].nunique()),'sample_rows':int(good.sum()),'minimum':float(valid.min()),'maximum':float(valid.max())}
 rec={'dataset':g,'all_expression_samples':len(m),'LM_expression_samples':len(z),'LM_patient_labels':z[idcol].nunique(),'dedup_key':idcol,'CGNL1_finite_LM_samples':finite,'survival_fields_in_GEO_matrix':clinical,'count_status':'verified','scope':'CRC liver metastases; strict colon-only subset not established in current matrix metadata'};bulk.append(rec)
 z.to_csv(A/'bulk'/(g+'_LM_only_metadata_CGNL1.tsv'),sep='\t',index=False)
(A/'bulk'/'BULK_LM_SCOPE_VERIFICATION.json').write_text(json.dumps(bulk,indent=2))
for g in ['GSE245552','GSE178318']:
 r=json.loads((A/(g+'_COVERAGE.json')).read_text());assert r['status']=='verified_counts_and_coverage'
assert json.loads((A/'GSE178318_ENDO_EXIT.json').read_text())['exit_code']==0
assert json.loads((A/'NEW_ENDO_EXIT.json').read_text())['exit_code']==0
D=Path(r'C:\Users\31937\Downloads\CRC肝转移_数据可用性核实_20260924');D.mkdir(exist_ok=True)
files=[]
for p in A.rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(A)
 if p.suffix in ['.md','.json','.py'] or (p.suffix=='.tsv' and not p.name.startswith('GPL')) or p.name.endswith('_metadata_excerpt.txt'):
  files.append((p,rel))
for p,rel in files:
 dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
manifest=[{'path':str(rel).replace('\\','/'),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p,rel in files]
(A/'AUDIT_DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(A/'AUDIT_DELIVERY_MANIFEST.json',D/'AUDIT_DELIVERY_MANIFEST.json')
archive=A/'cohort_verification_evidence_20260924.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as zf:
 for p,rel in files:zf.write(p,str(rel))
 zf.write(A/'AUDIT_DELIVERY_MANIFEST.json','AUDIT_DELIVERY_MANIFEST.json')
with zipfile.ZipFile(archive) as zf:assert zf.testzip() is None
shutil.copy2(archive,D/archive.name)
result={'audit_status':'complete','association_tests_run':False,'bulk_LM_patients':[171,167],'qualified_endo_total':14,'qualified_endo_external':9,'qualified_state_total':6,'qualified_state_external':3,'download_report':str(D/'数据可用性核实报告_20260924.md'),'evidence_zip_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'new_count_jobs_exit_code':0}
(A/'AUDIT_COMPLETE.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(A/'AUDIT_COMPLETE.json',D/'AUDIT_COMPLETE.json');print(json.dumps(result,ensure_ascii=False,indent=2));print(cs.to_string(index=False));print(json.dumps(bulk,ensure_ascii=False,indent=2))
