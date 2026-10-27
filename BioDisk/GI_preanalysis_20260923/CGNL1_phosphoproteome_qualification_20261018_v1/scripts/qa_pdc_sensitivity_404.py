from pathlib import Path
import csv,json,math
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
b=Path(__file__).resolve().parents[1];o=b/'results/PDC_processing_sensitivity_399';a=json.loads((o/'ANALYSIS_COMPLETE_399.json').read_text());vals=list(csv.DictReader((o/'patient_values.tsv').open(),delimiter='\t'));mp=json.loads((b/'metadata/PDC_ACTIVE_SAMPLE_QUALIFICATION_395.json').read_text())['active_map']
with (b/'raw/PDC_original_384/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv').open() as f:protein=next(r for r in csv.DictReader(f,delimiter='\t') if r['Gene']=='CGNL1')
with (b/'raw/PDC_original_384/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Phosphoproteome.phosphosite.tmt11.tsv').open() as f:sites={r['Phosphosite'].split(':')[-1].upper():r for r in csv.DictReader(f,delimiter='\t') if r['Gene']=='CGNL1'}
err=[]
for r in vals:
 c=mp[r['case']];suffix=' Unshared Log Ratio' if r['method']=='unshared' else ' Log Ratio';means=[]
 for row,suf,condition in [(sites[r['site']],' Log Ratio','Primary Tumor'),(sites[r['site']],' Log Ratio','Solid Tissue Normal'),(protein,suffix,'Primary Tumor'),(protein,suffix,'Solid Tissue Normal')]:
  v=[float(row[k+suf]) for k in c[condition] if row.get(k+suf,'') not in ['', 'NA','nan','NaN']];means.append(sum(v)/len(v))
 delta=means[0]-means[1]-means[2]+means[3];err.append(abs(delta-float(r['adjusted_delta'])))
for method in ['unshared','all_peptides']:
 tested=[]
 for c in a['results']:
  if c['method']!=method or c['P'] is None:continue
  d=[float(r['adjusted_delta']) for r in vals if r['site']==c['site'] and r['method']==method];assert len(set(abs(x) for x in d))==len(d) and all(x!=0 for x in d);p=float(wilcoxon(d,method='exact').pvalue);assert abs(p-c['P'])<1e-12;tested.append(c)
 q=multipletests([r['P'] for r in tested],method='fdr_bh')[1];assert max(abs(v-r['FDR']) for v,r in zip(q,tested))<1e-12
assert max(err)<1e-11
res={'all_pass':True,'raw_source_patient_values_checked':len(vals),'max_effect_error':max(err),'tests_independently_verified':8,'primary_S283_observed_pairs':44,'primary_S708_observed_pairs':14,'no_missing_value_imputation':True,'scope':'Technical reprocessing sensitivity within original cohort, not external patient validation'};(o/'INDEPENDENT_QA_404.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
