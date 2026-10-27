from pathlib import Path
import json,re,math,statistics,datetime
import numpy as np
from scipy.stats import wilcoxon
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=B/'results/HCC_exploratory_572';D=B/'results/core_qualification_568';norm=json.loads((O/'NORMALIZATION572.json').read_text());expected=json.loads((O/'patient_site_values572.json').read_text());result=json.loads((O/'RESULT572.json').read_text());pr=json.loads((D/'HCC_protein.json').read_text());sr=json.loads((D/'HCC_phospho.json').read_text());ps=next(s for s in pr['sheets'] if s['name']=='1. 6512 proteins');ss=next(s for s in sr['sheets'] if s['name']=='1. phosphoproteomic')
def keyed(s,row):
 headers={re.sub(r'\d+$','',k):v for k,v in s['headers'][1]['cells'].items()};return {headers[re.sub(r'\d+$','',k)]:v for k,v in row['cells'].items() if isinstance(headers.get(re.sub(r'\d+$','',k)),str) and re.fullmatch(r'\d+[TN]',headers[re.sub(r'\d+$','',k)])}
protein=keyed(ps,ps['target_rows'][0]);peptide_by_row={r['row']:(list(r['cells'].values())[3],keyed(ss,r)) for r in ss['target_rows']};errors=[];maxerr=0.
for e in expected:
 patient=e['patient'];T=patient+'T';N=patient+'N';dprotein=math.log2(protein[T]/protein[N])-(norm['protein_medians'][T]-norm['protein_medians'][N]);d={}
 for rn in e['complete_peptide_rows']:
  sequence,values=peptide_by_row[rn];v=math.log2(values[T]/values[N])-(norm['phospho_medians'][T]-norm['phospho_medians'][N]);d.setdefault(sequence,[]).append(v)
 actual=statistics.median(statistics.median(v) for v in d.values())-dprotein;error=abs(actual-e['adjusted_delta']);maxerr=max(maxerr,error)
 if error>1e-10:errors.append(dict(patient=patient,site=e['site'],delta_error=error))
checks=[]
for site in result['results']:
 if site['P'] is None:continue
 d=np.array([e['adjusted_delta'] for e in expected if e['site']==site['site']]);ties=len(set(abs(d)))!=len(d);zeros=int((d==0).sum());p=float(wilcoxon(d,method='exact',alternative='two-sided').pvalue);ok=abs(p-site['P'])<1e-10 and not ties and zeros==0;checks.append(dict(site=site['site'],scipy_exact_P=p,custom_exact_P=site['P'],ties=ties,zeros=zeros,pass_check=ok))
 if not ok:errors.append(dict(site=site['site'],test='exact signed-rank discrepancy/tie handling needs review'))
s283=np.array([e['adjusted_delta'] for e in expected if e['site']=='S283']);loo=[float(np.delete(s283,i).mean()) for i in range(len(s283))]
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_pass=not errors,patient_site_rows_checked=len(expected),max_effect_error=maxerr,test_checks=checks,errors=errors,S283_leave_one_patient_out_mean_range=[min(loo),max(loo)],scope='Independent scalar reconstruction from separately extracted source row dictionaries vs dense-matrix calculation; scipy exact test vs custom integer subset sums. Normalization offsets reused, not independently recomputed here. Does not resolve imputation/localization provenance.')
(O/'INDEPENDENT_QA574.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
