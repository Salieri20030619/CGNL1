from pathlib import Path
import csv,json,hashlib,gzip,datetime,collections,math
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');F=O/'C04_full_fraction_qualification';P=O/'C04_search_pilot';A=F/'joint_identification'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==z
spec=read(F/'JOINT_IDENTIFICATION_SPEC_BEFORE_RESULTS.json');inputs=[];notapp=[];missing=[]
for name in spec['files']:
 stem=Path(name).stem
 if name!='phospho_fraction_1.raw' and not (F/'searches'/stem/'FILE_SEARCH_COMPLETE.json').exists():missing.append(name);continue
 for mode,pilot,prefix in [('HCD','run_v2_all_spectra','pilot_hcd'),('CID','run_v3_CID','pilot_cid')]:
  w=P/pilot if name=='phospho_fraction_1.raw' else F/'searches'/stem/mode;prefix=prefix if name=='phospho_fraction_1.raw' else 'psms'
  if (w/'NOT_APPLICABLE.json').exists():notapp.append(read(w/'NOT_APPLICABLE.json'));continue
  if not (w/'SEARCH_EXIT.json').exists() or read(w/'SEARCH_EXIT.json')['exit']!=0:missing.append(name+':'+mode);continue
  assert sha(w/'comet.params')==spec[mode+'_parameter_template_sha256'];f=w/(prefix+'.txt');assert f.exists();inputs.append(dict(name=name,mode=mode,family=('pY_' if name=='phospho_pY.raw' else 'ST_')+mode,path=str(f),sha256=sha(f)))
if missing:raise SystemExit('Not ready;no pooledq computed. Missing: '+','.join(missing))
A.mkdir(exist_ok=True);save(A/'INPUT_MANIFEST.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputs=inputs,not_applicable=notapp,complete_frozen_file_scope=True))
csv.field_size_limit(20000000);rows=[];seen=set()
for item in inputs:
 with Path(item['path']).open(encoding='utf-8') as h:
  version=h.readline();assert version.startswith('CometVersion 2026.03')
  for x in csv.DictReader(h,delimiter='\t'):
   assert int(x['num'])==1;k=(item['name'],item['mode'],x['scan'],x['charge']);assert k not in seen;seen.add(k)
   x.update(file=item['name'],mode=item['mode'],family=item['family'],e_value=float(x['e-value']),is_decoy='DECOY_' in x['protein'],is_contaminant='CON__' in x['protein'],is_phospho='79.966' in x['modifications']);assert math.isfinite(x['e_value']) and x['e_value']>=0;rows.append(x)
def qval(rr,key):
 groups=collections.defaultdict(list)
 for x in rr:groups[x['e_value']].append(x)
 t=d=0;ordered=[]
 for score,g in sorted(groups.items()):
  d+=sum(x['is_decoy'] for x in g);t+=sum(not x['is_decoy'] for x in g);ordered.append((g,min(1.,(d+1)/max(t,1))))
 q=1.
 for g,fdr in reversed(ordered):
  q=min(q,fdr)
  for x in g:x[key]=q
families={};cg=[]
for family in sorted(set(x['family'] for x in rows)):
 rr=[x for x in rows if x['family']==family];pp=[x for x in rr if x['is_phospho']];qval(rr,'q_global');qval(pp,'q_phospho')
 passed=[x for x in pp if not x['is_decoy'] and x['q_global']<=.01 and x['q_phospho']<=.01]
 families[family]=dict(PSMs=len(rr),targets=sum(not x['is_decoy'] for x in rr),decoys=sum(x['is_decoy'] for x in rr),phospho_PSMs=len(pp),dual_q_pass_phospho_PSMs=len(passed))
for x in rows:
 if not x['is_decoy'] and any(a in x['protein'] for a in ['|Q6AW69','|D3Z3H4|']):cg.append(x)
fields=['file','mode','family','scan','charge','e_value','is_decoy','is_phospho','is_contaminant','q_global','q_phospho','plain_peptide','modified_peptide','protein','ascorepro','ascore_sitescores']
with gzip.open(A/'PSM_QVALUES.tsv.gz','wt',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
save(A/'CGNL1_ALL_CANDIDATES.json',cg)
save(A/'RESULTS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),families=families,total_PSMs=len(rows),CGNL1_candidates=len(cg),R_QA_pending=True,PSM_q_not_site_or_protein_FDR=True,biological_n=8,new_biological_tests=0,whole_project_complete=False));print(json.dumps(families))
