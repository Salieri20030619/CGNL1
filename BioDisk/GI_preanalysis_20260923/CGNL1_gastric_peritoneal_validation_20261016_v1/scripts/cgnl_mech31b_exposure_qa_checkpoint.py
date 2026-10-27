from pathlib import Path
import gzip,csv,math,json,statistics,datetime
import numpy as np
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2'
with gzip.open(O/'raw/GSE223165_readcount_genename.txt.gz','rt',encoding='utf-8') as f:
 reader=csv.reader(f,delimiter='\t');h=next(reader);cols=h[1:42];totals=[0]*41;target=None;length=None
 for r in reader:
  assert len(r)==51
  vals=[int(v) for v in r[1:42]];totals=[a+b for a,b in zip(totals,vals)]
  if r[0]=='ENSG00000128849':assert target is None;target=vals;length=int(r[h.index('gene_length')]);assert r[h.index('gene_name')]=='CGNL1'
assert target is not None
expected={c:(n/t*1e6,n/t*1e9/length) for c,n,t in zip(cols,target,totals)};err=0.;tests=[]
for mode in ['CPM','reconstructedFPKM']:
 with (O/'results'/f'GSE223165_{mode}_target_samples.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
 kind='CPM' if mode=='CPM' else 'FPKM';ix=0 if mode=='CPM' else 1
 for r in rows:err=max(err,abs(float(r[kind])-expected[r['matrix_column']][ix]),abs(float(r['log2'+kind+'1'])-math.log2(expected[r['matrix_column']][ix]+1)))
 with (O/'results'/f'GSE223165_{mode}_target_donor_deltas.tsv').open(encoding='utf-8') as f:ds=list(csv.DictReader(f,delimiter='\t'))
 with (O/'results'/f'GSE223165_{mode}_target_contrasts.tsv').open(encoding='utf-8') as f:rs=list(csv.DictReader(f,delimiter='\t'))
 for r in rs:
  d=[float(t[r['contrast']]) for t in ds if t[r['contrast']]];pv=float(stats.wilcoxon(d,method='exact').pvalue);assert abs(pv-float(r['primary_signedrank_P']))<1e-12;tests.append({'mode':mode,'contrast':r['contrast'],'n':len(d),'P_verified':pv})
 ps=[float(r['primary_signedrank_P']) for r in rs];order=sorted(range(2),key=lambda i:ps[i]);expectedq=[min(1,ps[order[0]]*2,ps[order[1]]),min(1,ps[order[1]])]
 for j,i in enumerate(order):assert abs(expectedq[j]-float(rs[i]['primary_FDR_two_target_contrasts']))<1e-12
assert err<1e-10
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'independent_scalar_counts_normalization_max_error':err,'independent_scipy_signedrank_checks':tests,'BH_independently_checked':True,'author_group_FPKM_reconstruction_verified':json.loads((O/'results/FPKM_RECONSTRUCTION_QA.json').read_text())['all_pass']}
(O/'results/EXPOSURE_NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
p=O/'STATUS.json';s=json.loads(p.read_text(encoding='utf-8'));s.update(utc=qa['utc'],last_goal_turn_classification='progress: within-mesothelial exposure response measured and independently verified; input corruption isolated; author FPKM reconstructed and allgroupvalidated')
s['exposure_result']={'GSE223165':'7/7 mesothelial donors show lower CGNL1 after all3 ovarian ascites exposures; donor-average FPKMlog1 difference -1.129684, signedrankP.015625,FDR.03125;10CM donors difference-.267888,P.105469. Exploratory ascites host response, not gastric-specific mechanism.','input_integrity':'failed compressed sampleFPKM preserved; intact annotatedcounts and author groupFPKM used; reconstructed58735genes x5group means allvalidated maxerror3.27e-9','qualification':'41samples are10mesothelial donors, with7donors exposed to3shared ascites; not41 or21 independent patients'}
s['qualified_new_datasets']['HRA002712']['access']='Public supplemental package obtained/CRCvalidated; source workbooks are figure data and signatures, not full transcriptome; no CGNL1 identifier hits; raw remains controlled; no access request sent'
s['next_actions']=['Assess whether ascites-related CGNL1 decline tracks overall mesothelial identity loss or a distinguishable state using explicitly frozen posthoc comparator panels','Cross-check gastric tissue source with GSE308231 once live download completes; current transfer must not be duplicated solely for slowness','Inspect FTP retry completion and compare intact FPKM if recovered; no repeated same-failure loops','Qualify other already discovered disease or stimulus datasets, especially direct gastric material','Integrate response direction with protein localization and functional rescue evidence without implying gastric causation'];p.write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa))
