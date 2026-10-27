from pathlib import Path
import json,datetime,hashlib,shutil
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p):return json.loads(p.read_text())
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8');t.replace(p)
results={};qa={}
for cohort,folder,q in [('H01','H01_582','INDEPENDENT_QA583.json'),('B03','B03_recovered_590','INDEPENDENT_QA591.json'),('H02','H02_recovered_590','INDEPENDENT_QA591.json')]:
 o=B/'results'/folder;results[cohort]=read(o/'RESULT.json');qa[cohort]=read(o/q);assert qa[cohort]['all_pass']
members=read(B/'metadata/resume_581/MEMBERS_COMPLETE585.json');assert members['exit_code']==0 and members['all_required_members_verified']
glob=read(B/'metadata/H02_GLOBAL_EXTRACTION589.json');assert glob['exit_code']==0 and glob['all_entries_read_without_libarchive_error']
for r in members['files']:
 p=Path(r['file']);assert p.stat().st_size==r['bytes'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
resolution=dict(utc=now,issues=[
 dict(id='H01-Q01',status='resolved_for_matrix_eligibility',finding='160 clinical pairs; case109 absent from both omics;159 matrix-eligible donors; exclusion saved'),
 dict(id='H01-Q02',status='retained_scientific_limit',finding='S283 changes sign without extra centering; do not call stable replication; source log base not independently resolved'),
 dict(id='H02-Q01',status='resolved',finding='Independent global-proteome archive PXD006512 recovered and validated; phosphoenriched proteinGroups not used as total denominator'),
 dict(id='H02-Q02',status='resolved_for_prespecified_subset',finding='Author quantitative-analysis crosswalk and BCLC A/0 inclusion yield53 mapped early-stage donors; other raw samples explicitly excluded;38 S283 complete pairs'),
 dict(id='SOURCE-Q01',status='unresolved',finding='G03/E02/B02/C02 public processed matrices remain unavailable or unqualified; prior access/source evidence retained; no outside contact')])
save(B/'metadata/ISSUE_RESOLUTION592.json',resolution)
gate=dict(utc=now,whole_project_complete=False,decision='Retain PDAC discovery and all prior negative validation; do not promote to a validated pan-digestive-cancer or cancer-specific mechanism',new_results=results,numerical_QA=qa,interpretation=dict(H01='S283 negative in centered primary but positive without extra centering; normalization-sensitive',B03='S283 negative primary; same direction but P>=.05 without extra centering. S486 negative in both analyses',H02='Early-stage clinically mapped raw-source subset: protein-adjusted S283 and S297 negative in primary and uncentered sensitivity; opposite to prior PDAC discovery direction; does not establish cancer-specificity or causality'),source_qualification=str(B/'results/source_qualification_586/QUALIFICATION.json'),prior_gate=str(B/'results/PHOSPHOSITE_MAINLINE_GATE_421.json'))
save(B/'results/PHOSPHOSITE_CONTINUATION_GATE592.json',gate)
plan=read(B/'metadata/bulk_576/PLAN.json');old={x['id']:x for x in read(B/'metadata/bulk_576/EXIT.json')['results']}
registry=[]
for job in plan['jobs']:
 cid=job['id'];status=old[cid]['stage'];detail='Previous source/analysis retained'
 if cid in results:status='paired_analysis_and_numerical_QA_complete';detail='Exploratory; retain cohort-specific source and sensitivity limits'
 elif cid in ['B01','C03']:status='download_verified_endpoint_qualification_complete';detail='Available CGNL1 aggregate matrices do not establish a site-specific paired endpoint'
 elif cid in ['G03','E02','B02','C02']:status='public_processed_source_followup_required';detail='See SOURCE_SEARCH579.json; no duplicate blind retries'
 elif cid=='B04':status='descriptive_analysis_retained';detail='Four pairs; total protein unavailable; no inferential P'
 elif cid in ['E03','H04','G04','C04','P05','P06','C05']:status='download_and_source_qualification_complete';detail='See source_qualification_586; not all sources support a patient-level site endpoint'
 registry.append(dict(id=cid,cohort=job.get('cohort'),status=status,detail=detail))
assert len(registry)==26
checkpoint=dict(utc=now,stage='resumed_analysis_batch_completed_and_qualified',whole_project_complete=False,new_completed_paired_analyses=['H01','B03','H02'],all_three_numerical_QA_pass=True,registry=registry,download_recovery=dict(H02='Required original phospho ZIP members and independent total-proteome RAR recovered, validated and analyzed',original_full_zip_complete=False,original_full_zip_parts_preserved=True),remaining_sources=['G03','E02','B02','C02'],scientific_gate=str(B/'results/PHOSPHOSITE_CONTINUATION_GATE592.json'),issue_resolution=str(B/'metadata/ISSUE_RESOLUTION592.json'),current_recovery_worker_running=False,legacy_monitor='Old monitor remains separate; its registered failed jobs are historical, not the completion evidence for this batch',next_actions=['Resolve public processed source availability for G03/E02/B02/C02 only when a new access route or source is identified','Retain normalization-sensitive H01 and source boundaries in cross-cohort synthesis','Do not infer cancer specificity from opposite direction or differing cohort P values'])
save(B/'metadata/PROJECT_CHECKPOINT_592.json',checkpoint)
statepath=B.parent/'CGNL1_CURRENT_STATE_20261017.json';backup=B/'metadata/CGNL1_CURRENT_STATE_BEFORE592.json'
assert not backup.exists();shutil.copyfile(statepath,backup);state=read(statepath)
state.update(utc=now,stage=checkpoint['stage'],whole_project_complete=False,active_root=str(B),latest_checkpoint=str(B/'metadata/PROJECT_CHECKPOINT_592.json'),latest_scientific_gate=checkpoint['scientific_gate'],next_native_completion=None,pending_launches=[],current_recovery_workers=[],remaining_public_source_followup=checkpoint['remaining_sources'])
state.pop('runtime_goal_state',None);state.pop('runtime_goal_state_verified_utc',None)
save(statepath,state)
artifacts=[B/'metadata/PROJECT_CHECKPOINT_592.json',B/'metadata/ISSUE_RESOLUTION592.json',B/'results/PHOSPHOSITE_CONTINUATION_GATE592.json',statepath]
save(B/'metadata/RESUME_ACCEPTANCE592.json',dict(utc=now,recovery_batch_complete=True,whole_project_complete=False,checks=['H01/B03/H02 numerical QA passed','Phospho member sizes CRC and SHA256 verified','Global RAR expected size SHA256 and complete libarchive read verified','All26 source statuses reconciled','Prior state backed up; prior validated scientific results retained'],artifacts=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in artifacts]))
print(json.dumps(dict(checkpoint=str(B/'metadata/PROJECT_CHECKPOINT_592.json'),recovery_batch_complete=True,whole_project_complete=False)))
