from pathlib import Path
import json,csv,decimal
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
b=Path(__file__).resolve().parents[1];o=b/'results/PDAC_phosphosite_357';a=json.loads((o/'ANALYSIS_COMPLETE_357.json').read_text());tables={}
for p in (b/'raw').glob('PDAC_*.txt'):
 if not any(k in p.name for k in ['proteomics_gene','phospho_site']):continue
 with p.open() as f:
  rd=csv.reader(f,delimiter='\t');header=next(rd);assert len(header)==len(set(header));rows={}
  for r in rd:
   if r[0].split('|')[0].split('.')[0]=='ENSG00000128849':rows[r[0]]=dict(zip(header[1:],r[1:]))
 tables[('site' if 'phospho' in p.name else 'protein','N' if '_Normal.' in p.name else 'T')]=rows
pt=tables['protein','T']['ENSG00000128849.10'];pn=tables['protein','N']['ENSG00000128849.10'];records=list(csv.DictReader((o/'patient_site_values.tsv').open(),delimiter='\t'));errs=[];checks=[];family=[]
for c in a['contrasts']:
 rows=[r for r in records if r['site']==c['site']];ds=[]
 for r in rows:
  k=r['case'];v=[tables['site','T'][c['site']][k],tables['site','N'][c['site']][k],pt[k],pn[k]];q=[decimal.Decimal(z) for z in v];d=float(q[0]-q[1]-q[2]+q[3]);ds.append(d);errs.append(abs(d-float(r['adjusted_delta'])))
 assert len(ds)==c['n_pairs']
 if len(ds)>=10:
  assert len(set(abs(v) for v in ds))==len(ds) and all(v!=0 for v in ds)
  p=float(wilcoxon(ds,alternative='two-sided',method='exact').pvalue);assert abs(p-c['P'])<1e-12;family.append(c);checks.append({'site':c['site'].split('|')[2],'n':len(ds),'P':p})
q=multipletests([c['P'] for c in family],method='fdr_bh')[1];assert max(abs(v-c['FDR']) for v,c in zip(q,family))<1e-12;assert max(errs)<1e-11
res={'all_pass':True,'paired_cases':len(a['cases']),'source_decimal_patient_site_values_checked':len(records),'site_tests_independently_checked':len(checks),'max_delta_roundoff':max(errs),'checks':checks,'scope':'Numerical/source-row QA; source imputation provenance and independent cohort replication remain pending'};(o/'INDEPENDENT_QA_360.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in res.items() if k!='checks'}))
