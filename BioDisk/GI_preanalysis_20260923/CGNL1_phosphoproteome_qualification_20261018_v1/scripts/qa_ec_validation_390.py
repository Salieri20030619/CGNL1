from pathlib import Path
import json,csv,math,openpyxl
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
b=Path(__file__).resolve().parents[1];o=b/'results/EC_phosphosite_validation_388';a=json.loads((o/'ANALYSIS_COMPLETE_388.json').read_text());records=list(csv.DictReader((o/'patient_site_values.tsv').open(),delimiter='\t'));w=openpyxl.load_workbook(b/'raw/41467_2021_25202_MOESM3_ESM.xlsx',read_only=True,data_only=True);source={}
for sn,col,offset in [('b. 9,300 proteins',1,2),('c. 61,471 phosphosites',2,4)]:
 it=w[sn].iter_rows(values_only=True);h=next(it)
 for r in it:
  if r[col]=='CGNL1':source['protein' if offset==2 else r[3]]=dict(zip(h[offset:],r[offset:]))
 if offset==2:assert 'protein' in source
errors=[];checks=[]
for r in records:
 site=r['feature'].split('/')[0];case=r['patient'];p=source['protein'];v=source[site];vt=v[case+'_T'];vn=v[case+'_N'];pt=p[case+'T'];pn=p[case+'N'];assert vt>0 and vn>0 and vt!=0.01 and vn!=0.01;d=(math.log(vt)-math.log(vn))/math.log(2)-pt+pn;errors.append(abs(d-float(r['adjusted_delta'])))
for c in a['results']:
 label='/'.join(c['feature_labels']);ds=[float(r['adjusted_delta']) for r in records if r['feature']==label];assert len(ds)==c['n_pairs']
 for alias in c['feature_labels']:assert source[alias]==source[c['feature_labels'][0]]
 if c['P'] is not None:
  assert len(set(abs(v) for v in ds))==len(ds) and all(v!=0 for v in ds);p=float(wilcoxon(ds,method='exact').pvalue);assert abs(p-c['P'])<1e-12;checks.append(c)
qs=multipletests([c['P'] for c in checks],method='fdr_bh')[1];assert max(abs(q-c['FDR']) for q,c in zip(qs,checks))<1e-12 and max(errors)<1e-11
res={'all_pass':True,'patient_feature_values_checked_from_original_XLSX':len(records),'unique_feature_tests':len(checks),'max_effect_error':max(errors),'site_alias_vectors_verified':True,'scope':'Observed-pair numerical verification and original-source rows; cross-disease validation not PDAC-internal replication'};(o/'INDEPENDENT_QA_390.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
