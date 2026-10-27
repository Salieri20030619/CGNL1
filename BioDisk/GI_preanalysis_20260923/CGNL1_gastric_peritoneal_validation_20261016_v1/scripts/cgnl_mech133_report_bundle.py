from pathlib import Path
import json,zipfile,hashlib,datetime
O=Path.cwd()/'mechanism_investigation_v2';dest=O/'report/review_update_v2';dest.mkdir(exist_ok=True);files=[]
for prefix,root,names in [('baseline',O.parent/'results',['ORIGINAL31_RESULT.json','BULK_RESULT.json','original31_models.tsv','bulk_models.tsv','identity_source_contrasts.tsv','mesothelial_state_models.tsv','bulk_patient_deltas.tsv']),('current',O/'results',['MECHANISM_EVIDENCE_LEDGER_v8.json','GSE237876_complete_route_summary.tsv','GSE237876_paired_deltas.tsv','GSE237876_sample_values.tsv','ROUTE_NUMERICAL_QA_100.json','GSE314812_GC05_CPM_sensitivity.tsv','GC05_NUMERICAL_QA_110.json','GSE223165_reconstructedFPKM_target_contrasts.tsv','GSE223165_reconstructedFPKM_target_donor_deltas.tsv','GSE223165_identity_comparator_tests.tsv','AUTHOR_FPKM_RECONSTRUCTION_CROSSCHECK.json','GSE53276_contrasts.tsv','HPA_GASTRIC_PROTEIN_SUMMARY.json'])]:
 for name in names:
  p=root/name;assert p.is_file();files.append((prefix+'/'+name,p))
for name in ['CLINICAL_QUALIFICATION_109.json','GSE308231_SPECIES_AND_METADATA_106.json','COHORT_PROVENANCE_102.json','CLINICAL_PROVENANCE_103.json','SPATIAL_FEASIBILITY_124.json','PIPS_SPATIAL_FINAL_QUALIFICATION_131.json','LOCAL_ACCELERATION_127.json']:
 p=O/'metadata'/name;assert p.is_file();files.append(('qualification/'+name,p))
files.extend([('specification/MECHANISM_MAINLINE_DECISION_CRITERIA_118.json',O/'MECHANISM_MAINLINE_DECISION_CRITERIA_118.json'),('specification/DRA_TARGET_STATISTICAL_SPECIFICATION.json',O/'DRA_TARGET_STATISTICAL_SPECIFICATION.json')])
manifest={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'review update, not final mechanism report','files':{name:{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'source':str(p)} for name,p in files}}
p=dest/'evidence_bundle.zip'
with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED) as z:
 for name,source in files:z.write(source,name)
 z.writestr('MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2).encode('utf-8'))
with zipfile.ZipFile(p) as z:assert z.testzip() is None
print(json.dumps({'zip':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'entries':len(files)+1,'bytes':p.stat().st_size}))
