
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,csv,gzip,math,datetime,subprocess
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
p=Path.cwd();old=p.parent/'CGNL1_gastric_peritoneal_validation_20261016_v1'
launch=json.loads((p/'metadata/SPHEROID_ANALYSIS_LAUNCH_233.json').read_text())
if not Path(launch['receipt']).exists():subprocess.run(['/usr/bin/python3',str(old/'scripts/cgnl_wait_pidfd_system.py'),str(launch['wrapper_pid'])],check=True)
assert Path(launch['receipt']).exists() and json.loads(Path(launch['receipt']).read_text())['exit_code']==0
out=p/'results/spheroid_233';complete=json.loads((out/'ANALYSIS_COMPLETE_233.json').read_text())
tasks=json.loads((p/'metadata/DOWNLOAD_PLAN_231.json').read_text())['tasks'];shared=set(json.loads((out/'COMMON_GENE_UNIVERSE.json').read_text()))
saved=pd.read_csv(out/'sample_target_summary.tsv',sep=chr(9)).set_index('accession');checks=[];rawtotals={}
for task in tasks:
 acc=task['accession'];file=p/'raw'/task['name'];totals={};special={}
 with gzip.open(file,'rt',encoding='utf-8',newline='') as h:
  reader=csv.reader(h,delimiter=chr(9));header=next(reader);first=next(reader)
  barcodes=header[1:] if len(first)==len(header) else header
  import itertools
  for row in itertools.chain([first],reader):
   assert len(row)==len(barcodes)+1 and row[0] not in totals
   values=[int(float(v)) for v in row[1:]];totals[row[0]]=sum(values)
   if row[0] in ['CGNL1','MKI67','TOP2A']:special[row[0]]=values
 with gzip.open(p/'processed/spheroid_233'/(acc+'_gene_totals.json.gz'),'rt',encoding='utf-8') as h:prior=json.load(h)
 assert totals==prior
 arrays=np.load(p/'processed/spheroid_233'/(acc+'_cell_counts.npz'),allow_pickle=False);assert arrays['cell_id'].tolist()==barcodes
 for g,v in special.items():assert arrays[g].tolist()==v
 r=saved.loc[acc];assert len(barcodes)==r.cells and sum(totals.values())==r.provided_UMI
 denominator=sum(totals[g] for g in shared);assert denominator==r.common_UMI
 present='CGNL1' in totals;assert present==bool(r.target_measured)
 if present:
  assert totals['CGNL1']==r.target_UMI and sum(v>0 for v in special['CGNL1'])==r.target_detected_cells
  assert abs(math.log2(1+totals['CGNL1']/denominator*1e6)-r.CGNL1_common_log2CPM1)<1e-12
  assert abs(math.log2(1+totals['CGNL1']/sum(totals.values())*1e6)-r.CGNL1_provided_log2CPM1)<1e-12
 rawtotals[acc]=totals;checks.append({'accession':acc,'gene_totals_exact':True,'cells':len(barcodes),'target_measured':present})
 print('INDEPENDENT_RAW_QA_PASS',acc,flush=True)
assert set.intersection(*(set(v) for v in rawtotals.values()))==shared
if complete['all_target_measured']:
 wide=saved.reset_index().pivot(index='patient',columns='condition',values='CGNL1_common_log2CPM1');delta=wide.spheroid-wide.parent
 primary=next(r for r in complete['contrasts'] if r['role']=='primary')
 assert len(delta)==3 and abs(float(delta.mean())-primary['mean_difference'])<1e-12
 if (delta!=0).all():assert abs(float(wilcoxon(delta,method='exact').pvalue)-primary['P'])<1e-12
else:assert complete['contrasts'][0]['status']=='not_evaluable_due_to_source_filtered_target_absence'
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'input_matrices_checked':6,'checks':checks,'method':'Independent Python integer accumulation of every raw matrix gene row; every pseudobulk gene total, barcodes and target/cell-cycle raw vectors checked against first implementation; normalization and primary exact test checked','scientific_interpretation_pending':True,'whole_project_complete':False}
(out/'INDEPENDENT_RAW_QA_234.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2),flush=True)
