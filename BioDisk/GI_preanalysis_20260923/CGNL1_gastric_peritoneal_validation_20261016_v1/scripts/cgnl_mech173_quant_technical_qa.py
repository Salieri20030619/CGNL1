from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,math,datetime,pandas as pd,numpy as np,hashlib
O=Path.cwd()/'mechanism_investigation_v2';state=json.loads((O/'metadata/DRA_QUANTIFICATION_STATUS.json').read_text(encoding='utf-8'));idx=json.loads((O/'metadata/INDEX_COMPLETE.json').read_text(encoding='utf-8'));out=[]
for run,rec in state['done'].items():
 assert rec['exit_code']==0;p=O/'processed/DRA014736_quant'/run;meta=json.loads((p/'aux_info/meta_info.json').read_text(encoding='utf-8'));assert meta['quant_errors']==[] and meta['index_seq_hash']==idx['info']['SeqHash'];x=pd.read_csv(p/'quant.sf',sep='\t');g=pd.read_csv(p/'gene_quant.tsv.gz',sep='\t');assert x.Name.is_unique and not g[['gene_id','gene_symbol']].duplicated().any();assert np.isfinite(x[['TPM','NumReads']].values).all() and np.isfinite(g[['TPM','NumReads']].values).all();totals={};errors={}
 for col in ['TPM','NumReads']:
  assert (x[col]>=0).all() and (g[col]>=0).all();tx=math.fsum(x[col]);gg=math.fsum(g[col]);errors[col]=abs(tx-gg);assert errors[col]<1e-5;totals[col]=tx
 assert abs(totals['TPM']-1e6)<1.;target=g[g.gene_id.str.split('.').str[0]=='ENSG00000128849'];assert len(target)==1 and target.gene_symbol.iloc[0]=='CGNL1';out.append({'run':run,'donor':rec['donor'],'condition':rec['condition'],'mapping_percent':meta['percent_mapped'],'library_types':meta['library_types'],'gene_transcript_sum_error':errors,'total_transcript_TPM':totals['TPM'],'estimated_reads_minus_num_mapped':totals['NumReads']-meta['num_mapped'],'target_unique_row_present':True,'gene_quant_sha256':hashlib.sha256((p/'gene_quant.tsv.gz').read_bytes()).hexdigest()})
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completed_runs_checked':len(out),'expected_runs':12,'complete':len(out)==12 and not state['fail'],'technical_checks_pass_for_completed_runs':True,'runs':out,'target_response_outcomes_not_evaluated':True};(O/'results/DRA_TECHNICAL_QA_173.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
