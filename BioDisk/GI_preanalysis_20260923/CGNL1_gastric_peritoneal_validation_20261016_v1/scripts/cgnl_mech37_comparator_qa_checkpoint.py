from pathlib import Path
import json,gzip,csv,math,statistics,datetime
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2';spec=json.loads((O/'IDENTITY_COMPARATOR_SPECIFICATION.json').read_text());panels=spec['panels'];wanted=set(sum(panels.values(),[])+['CGNL1']);values={};sums=[0]*41
with gzip.open(O/'raw/GSE223165_readcount_genename.txt.gz','rt',encoding='utf-8') as f:
 rd=csv.reader(f,delimiter='\t');h=next(rd);cols=h[1:42]
 for r in rd:
  v=[int(x) for x in r[1:42]];sums=[a+b for a,b in zip(sums,v)];g=r[h.index('gene_name')]
  if g in wanted:assert g not in values;values[g]=(v,int(r[h.index('gene_length')]))
normalized={g:[math.log2(n/t*1e9/l+1) for n,t in zip(v,sums)] for g,(v,l) in values.items()};maxerr=0
with (O/'results/GSE223165_identity_comparator_scores.tsv').open(encoding='utf-8') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  i=cols.index(r['matrix_column'])
  for panel,gs in panels.items():maxerr=max(maxerr,abs(float(r[panel])-statistics.mean(normalized[g][i] for g in gs)))
assert maxerr<1e-10
with (O/'results/GSE223165_identity_comparator_deltas.tsv').open(encoding='utf-8') as f:ds=list(csv.DictReader(f,delimiter='\t'))
with (O/'results/GSE223165_identity_comparator_tests.tsv').open(encoding='utf-8') as f:rs=list(csv.DictReader(f,delimiter='\t'))
for r in rs:
 d=[float(x['difference']) for x in ds if x['measure']==r['measure'] and x['exposure']==r['exposure']];assert abs(stats.wilcoxon(d,method='exact').pvalue-float(r['P']))<1e-12
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'source_count_to_all_panel_values_independently_recomputed':True,'max_error':maxerr,'signedrank_comparisons_verified':len(rs)};(O/'results/IDENTITY_COMPARATOR_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
ledger=json.loads((O/'results/MECHANISM_EVIDENCE_LEDGER_v1.json').read_text(encoding='utf-8'))
ledger.append({'id':'MEC-05','candidate':'Ascites-associated CGNL1 decline is a distinct mesothelial functional state beyond lineage attenuation','evidence':'Human ovarian ascites exposure:7/7donors targetlower, but target-minus5marker identity difference -.0288,P.6875; target-minus2specificmarkers difference+.0867,P1.0','interpretation':'Current RNA data do not show preferential target attenuation beyond lineage marker decline. Does not prove biological equivalence or absence of function.','positive_control':'CM in10donors increases mesenchymal score in10/10 and lowers junction score9/10, while target changeP.1055. Ascites increases ANGPTL4/STC1 in7/7 descriptively.','decision':'Do not promote CGNL1 as MMT driver on these data; prioritize direct gastric protein/source evidence','data':'GSE223165','analysis':'posthoc identity comparator, family sizes4 and6; donor-level units; ovarian exposure not gastric-specific'})
(O/'results/MECHANISM_EVIDENCE_LEDGER_v2.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2),encoding='utf-8')
p=O/'STATUS.json';s=json.loads(p.read_text(encoding='utf-8'));s.update(utc=qa['utc'],last_goal_turn_classification='progress: lineage comparator analysis and independent QA completed; specific-driver interpretation weakened; directGC mesothelial experimental data accession DRA014736 identified',mechanism_ledger='results/MECHANISM_EVIDENCE_LEDGER_v2.json');s['direct_gastric_followup']={'study':'Yonemura et al Cell Reports2024 DOI10.1016/j.celrep.2023.113613','reported_RNA_accession':'DRA014736','source':'https://kumadai.repo.nii.ac.jp/record/2000526/files/igaku_kou2464zenbun.pdf','status':'ENA run metadata request in progress; verify primary design before download','caveat':'MMT/TNC colonization axis already established in that paper, not novel as a generic pathway'};p.write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa))
