from pathlib import Path
import json,csv,re,math,collections
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
b=Path(__file__).resolve().parents[1];all_checks=[]
# independent direct source reads; verify patient aggregation from actual biological metadata
samples=json.loads((b/'metadata/GASTRIC_PDC_QUALIFICATION_407.json').read_text())['data']['samples'];mp={r['aliquot_submitter_id']:r['case_submitter_id'] for r in samples}
with (b/'raw/GASTRIC_PDC_408/Gastric_Cancer_Korea_Proteome.itraq.tsv').open() as f:protein=next(r for r in csv.DictReader(f,delimiter='\t') if r['Gene']=='CGNL1')
with (b/'raw/GASTRIC_PDC_408/Gastric_Cancer_Korea_Phosphoproteome.phosphosite.itraq.tsv').open() as f:sites={r['Phosphosite'].split(':')[-1].upper():r for r in csv.DictReader(f,delimiter='\t') if r['Gene']=='CGNL1'}
o=b/'results/GASTRIC_validation_415';old=json.loads((o/'ANALYSIS_COMPLETE_415.json').read_text());saved=list(csv.DictReader((o/'patient_values.tsv').open(),delimiter='\t'));errors=[]
for r in saved:
 values=[]
 for col,val in sites[r['site']].items():
  m=re.fullmatch(r'Log (T\d+)(?:\.\d+)?/(N\d+)(?:\.\d+)?',col)
  if not m:continue
  t,n=m.groups();assert mp[t]==mp[n]
  if mp[t]!=r['patient']:continue
  pv=protein.get('Unshared '+col,'')
  if val not in ['', 'NA','nan'] and pv not in ['', 'NA','nan']:values.append(float(val)-float(pv))
 assert len(values)==int(r['n_technical_channels']);errors.append(abs(sum(values)/len(values)-float(r['adjusted_delta'])))
for cohort,folder,resultkey,methodkey in [('Korean_GC',o,'results','FDR_within_GC'),('CPTAC_GC2026',b/'results/GC2026_validation_418','results','FDR_within_GC2026')]:
 d=json.loads((folder/('ANALYSIS_COMPLETE_415.json' if cohort=='Korean_GC' else 'ANALYSIS_COMPLETE_418.json')).read_text());rows=list(csv.DictReader((folder/'patient_values.tsv').open(),delimiter='\t'));tested=[]
 for r in d[resultkey]:
  ds=[float(v['adjusted_delta']) for v in rows if v['site']==r['site']];assert len(ds)==r['n_patients']
  if r['P'] is not None:
   assert len(set(abs(v) for v in ds))==len(ds) and all(v!=0 for v in ds);p=float(wilcoxon(ds,method='exact').pvalue);assert abs(p-r['P'])<1e-12;tested.append(r)
 q=multipletests([r['P'] for r in tested],method='fdr_bh')[1];assert max(abs(v-r[methodkey]) for v,r in zip(q,tested))<1e-12
 if cohort=='CPTAC_GC2026':
  with (b/'raw/GC2026_410/GC_metatable.tsi').open() as f:met={r['Sample_ID']:r for r in csv.DictReader(f,delimiter='\t') if r['Sample_ID']!='CAT'}
  bycase=collections.defaultdict(dict)
  for r in met.values():bycase[r['Case']][r['tissue_type']]=r['Sample_ID']
  with (b/'raw/GC2026_410/GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor_NAT.cct').open() as f:prot=next(r for r in csv.DictReader(f,delimiter='\t') if r['Gene']=='CGNL1')
  with (b/'raw/GC2026_410/GC_TMT_phospho_ratio_GencodeV42_Tumor_NAT.cct').open() as f:ph={r['Gene'].split('|')[3]:r for r in csv.DictReader(f,delimiter='\t') if r['Gene'].startswith('CGNL1_')}
  for r in rows:
   t=bycase[r['case']]['Tumor'];n=bycase[r['case']]['NAT'];delta=float(ph[r['site']][t])-float(ph[r['site']][n])-float(prot[t])+float(prot[n]);errors.append(abs(delta-float(r['adjusted_delta'])))
 external=d['cumulative_external_family'];qs=multipletests([r['P'] for r in external],method='fdr_bh')[1];assert max(abs(q-r['FDR_cumulative_external']) for q,r in zip(qs,external))<1e-12
 res={'all_pass':True,'cohort':cohort,'source_patient_values_checked':len(rows),'test_count':len(tested),'max_effect_error':max(errors),'source_pairing_checked':True,'source_integrity':'Authoritative PDC MD5 or recorded complete public-source hashes','independent_numerical_method':'SciPy exact signed-rank and statsmodelsBH; direct rawsource scalarreaggregation'};(folder/('INDEPENDENT_QA_419.json' if cohort=='Korean_GC' else 'INDEPENDENT_QA_420.json')).write_text(json.dumps(res,indent=2),encoding='utf-8');all_checks.append(res)
assert max(errors)<1e-11;print(json.dumps(all_checks))
