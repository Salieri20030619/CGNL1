from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,os
root=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
out=root/'spatial_mainline_validation_20261016_v1'
out.mkdir(exist_ok=False)
for name in ['baseline','scripts','logs','metadata','literature','results','figures','report','downloads']:(out/name).mkdir()
def write(name,obj): (out/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
record={
 'created_utc':datetime.now(timezone.utc).isoformat(),'host':os.uname().nodename,
 'user_authorization':'先记录当前状态，以这个基点自行进行验证，允许下载数据和服务器上进行任何方法，保留完整的记录，到核定最终主线报告让我审核，可行性和创新性都需要纳入考虑',
 'objective':'Evaluate and empirically validate a worthwhile spatial research question for CGNL1 in colorectal cancer liver metastasis, assessing feasibility and novelty before selecting a final mainline for user review.',
 'candidate_question':'Does CGNL1 mark an epithelial state associated with tumor-liver interface organization, beyond generic epithelial differentiation and local cell composition?',
 'candidate_not_presumed_true':True,'alternative_explanations':['tissue and lineage mixture','generic mature epithelial state','regional technical quality','patient-specific organization','spatial capture and sampling'],
 'prior_state':'SCP2079 paired validation complete; entire CGNL1 project incomplete. Prior broad RNA-based stop recommendation was too broad for spatial question. M242T technical pilot lacks CGNL1 measurement.',
 'known_results_before_new_design':['original23 actomyosin M1 negative P0.01119 FDR0.03358, external inconsistency','GSE294385 mural neighborhood signal attenuated in full adjustment','7-patient micro/macro CGNL1 contrast attenuated after state and mixture adjustment','SCP2079 n14 decrease in11 P0.01660 FDR0.04980; prediction sensitivity failure'],
 'analysis_policy':{'unit':'independent biological patient','nominal_gate':'target-effect raw P<0.05','FDR':'report within predefined comparable hypothesis families; not sole selection gate','small_n':'exploratory direction, effect, concordance; no inference without valid biological replication','discovery_status':'existing cohorts already explored; new estimands on these patients are exploratory, not untouched external validation','no_causal_claim_from_cross_section':'required'},
 'preserve_originals':True,'local_policy':'analysis data/code/intermediates remain on sgs10; local final reports and handoff only; temporary upload/orchestration helpers are separate',
 'deliverable':'full traceable mainline decision report for user review; not automatic manuscript/client release','external_messages_authorized':False,
 'status':'baseline_recorded'
}
write('BASELINE_AND_AUTHORIZATION.json',record)
files=['CGNL1_CURRENT_STATE_20261014.json','CGNL1_PROGRESS_20261014_SCP2079_PAIRED.json','V3_PROGRESS_20261007_FULLSLICE.json','epithelial_independence_20261012_v1/INDEPENDENCE_DECISION.json','reassessment_20261008_v1/results/R1_target_summary.tsv','feasibility_audit_20261008_v1/results/GSE294385_v3_patient_summary.tsv','epithelial_independence_20261012_v1/results/platform_adapted/C_gene_summary.tsv']
manifest=[]
for rel in files:
 p=root/rel; b=p.read_bytes();dest=out/'baseline'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b);manifest.append({'source':str(p),'snapshot':str(dest.relative_to(out)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
write('baseline/SOURCE_MANIFEST.json',manifest)
write('STATUS.json',{'stage':'baseline_recorded','utc':record['created_utc'],'mainline_decision':'pending'})
print(json.dumps({'directory':str(out),'snapshots':len(manifest),'baseline_recorded':True}))
