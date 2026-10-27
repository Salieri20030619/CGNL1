from pathlib import Path
import json,sys,hashlib,csv,gzip,collections,datetime,math
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1')
name=sys.argv[1];assert name in ['protein_set1_fraction_1.raw','protein_set2_fraction_1.raw'];w=O/'C04_totalprotein_pilot_20261026'/Path(name).stem/'search_CID'
r=json.loads((w/'SEARCH_EXIT.json').read_text());assert r['exit']==0
for item in r['files']:
 f=w/item['name'];assert f.stat().st_size==item['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==item['sha256']
assert not (w/'PSM_QA_SPEC.json').exists()
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=name+' pilot only',ordering='E-value ascending, ties grouped',q_estimator='reverse cumulative min of (decoys+1)/max(targets,1), capped1',q_threshold=.01,no_protein_FDR_claim=True,no_absence_inference=True,no_disease_tests=True)
(w/'PSM_QA_SPEC.json').write_text(json.dumps(spec,indent=2)+'\n',encoding='utf-8')
genes=set()
for line in (O/'C04_search_pilot/mouse_plus_contaminants.fasta').open():
 if line.startswith('>') and 'GN=Cgnl1 ' in line:genes.add(line[1:].split()[0])
rows=[];seen=set();csv.field_size_limit(20000000)
with (w/'psms.txt').open() as h:
 assert h.readline().startswith('CometVersion 2026.03')
 for x in csv.DictReader(h,delimiter='\t'):
  assert int(x['num'])==1;key=(x['scan'],x['charge']);assert key not in seen;seen.add(key);assert '79.966' not in x['modifications'];e=float(x['e-value']);assert math.isfinite(e) and e>=0
  rows.append(dict(scan=int(x['scan']),charge=int(x['charge']),e_value=e,is_decoy='DECOY_' in x['protein'],is_contaminant='CON__' in x['protein'],peptide=x['plain_peptide'],protein=x['protein']))
groups=collections.defaultdict(list)
for x in rows:groups[x['e_value']].append(x)
t=d=0;ordered=[]
for score,g in sorted(groups.items()):
 d+=sum(x['is_decoy'] for x in g);t+=sum(not x['is_decoy'] for x in g);ordered.append((g,min(1,(d+1)/max(t,1))))
q=1
for g,v in reversed(ordered):
 q=min(q,v)
 for x in g:x['q_pilot']=q
cg=[]
for x in rows:
 if not x['is_decoy'] and genes.intersection(x['protein'].split(',')):
  z=dict(x);z['all_reported_proteins_Cgnl1']=set(x['protein'].split(','))<=genes;cg.append(z)
with gzip.open(w/'PSM_QVALUES.tsv.gz','wt',encoding='utf-8',newline='') as h:
 writer=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
(w/'CGNL1_CANDIDATES.json').write_text(json.dumps(cg,indent=2)+'\n',encoding='utf-8')
passed=[x for x in rows if not x['is_decoy'] and not x['is_contaminant'] and x['q_pilot']<=.01]
r=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),PSMs=len(rows),targets=t,decoys=d,noncontaminant_target_PSM_q01=len(passed),unique_passing_peptides=len({x['peptide'] for x in passed}),CGNL1_candidate_PSMs=len(cg),CGNL1_pilot_q01=sum(x['q_pilot']<=.01 and not x['is_contaminant'] for x in cg),all_search_output_hashes_verified=True,independent_R_QA_pending=True,protein_FDR_not_established=True,full_cohort_not_run=True,whole_project_complete=False)
(w/'PSM_QA_RESULTS.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))
