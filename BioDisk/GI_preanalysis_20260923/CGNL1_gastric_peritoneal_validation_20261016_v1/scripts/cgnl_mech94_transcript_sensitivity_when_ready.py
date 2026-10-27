from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,datetime
import pandas as pd,numpy as np
O=Path.cwd()/'mechanism_investigation_v2';ann=json.loads((O/'metadata/TARGET_TRANSCRIPT_BIOTYPES_93.json').read_text());target=[r for r in ann['transcripts'] if r['gene_symbol']=='CGNL1'];coding=ann['CGNL1_coding_transcripts'];ids=[r['transcript_id'] for r in target];assert len(ids)==3 and len(coding)==1
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'outcomes_not_yet_observed':True,'primary_gene_analysis_unchanged':True,'purpose':'separateannotatedCDS-definedCGNL1transcriptfromtwoCDS-not-definedtranscriptsbeforeproteinmechanisminterpretation','coding':coding,'other':[t for t in ids if t not in coding],'annotation_wording':'protein_coding_CDS_not_defined doesnotprove noncoding,functionless,oruntranslated; retainexactlabel','measure':'individualtranscriptTPMandestimatedNumReads;codingTPM/allCGNL1TPMwhenall>0;log2codingTPM+1pairedchanges','inference':'descriptivesensitivityonly,n3donors,noadditionalgene-selectionPvalues','mapping_uncertainty':'sharedreadsmaylimitisoformidentifiability; noisoformswitchingclaimwithoutadequateassignment evidence'};(O/'DRA_TRANSCRIPT_SENSITIVITY_SPECIFICATION.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
done=O/'metadata/DRA_QUANTIFICATION_COMPLETE.json'
def live():
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    if any(a.endswith(b'/cgnl_mech46b_quantify_when_ready.py') for a in (p/'cmdline').read_bytes().split(b'\0')):return True
   except (FileNotFoundError,PermissionError):pass
 return False
while not done.exists() and live():time.sleep(900)
assert done.exists() and json.loads(done.read_text())['all_pass'],'Incompletequantification'
plan=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text());runs={t['run']:t for t in plan['tasks']};rows=[];summary=[]
for run,t in runs.items():
 a=pd.read_csv(O/'processed/DRA014736_quant'/run/'quant.sf',sep='\t').set_index('Name');assert all(i in a.index for i in ids);total=float(a.loc[ids,'TPM'].sum());ct=float(a.loc[coding,'TPM'].sum());summary.append({'run':run,'donor':t['donor'],'condition':t['condition'],'all_gene_TPM':total,'defined_CDS_TPM':ct,'defined_CDS_fraction':ct/total if total>0 else np.nan,'log2definedCDSTPM1':np.log2(ct+1)})
 for r in target:rows.append(dict(r,run=run,donor=t['donor'],condition=t['condition'],TPM=float(a.loc[r['transcript_id'],'TPM']),estimated_NumReads=float(a.loc[r['transcript_id'],'NumReads'])))
s=pd.DataFrame(summary);v=s.pivot(index='donor',columns='condition',values='log2definedCDSTPM1');d=pd.DataFrame({c+'_minus_control':v[c]-v.control for c in ['combined','TGF','IL1']})
pd.DataFrame(rows).to_csv(O/'results/DRA014736_transcript_values.tsv',sep='\t',index=False,encoding='utf-8');s.to_csv(O/'results/DRA014736_defined_CDS_summary.tsv',sep='\t',index=False,encoding='utf-8');d.to_csv(O/'results/DRA014736_defined_CDS_donor_deltas.tsv',sep='\t',encoding='utf-8');(O/'results/DRA_TRANSCRIPT_SENSITIVITY_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'donors':3,'libraries':12,'transcripts':3,'P_values_computed':False,'protein_effect_not_inferred':True},indent=2),encoding='utf-8');print(d.to_string(),flush=True)
