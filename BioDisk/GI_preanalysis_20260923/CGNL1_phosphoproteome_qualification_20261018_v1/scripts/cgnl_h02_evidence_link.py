from pathlib import Path
import json,re,hashlib,datetime,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'H02_localization_reaudit_20261024';W=D/'evidence_link';W.mkdir(exist_ok=True)
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==z
src=D/'H02_CGNL1_evidence_rows.json';assert sha(src)==read(D/'H02_target_evidence_receipt.json')['sha256_target']
save(W/'SPEC_BEFORE_COUNTS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(src),scope='68originalH02measurements,17patients,S283/S297',linkage='exactExperiment and reciprocal siteID/evidenceID links;probabilitysequence uniquely mapscanonicalCGNL1',direct_rule='TypeMSMSorMULTI-MSMS withMS/MSCount>0 andnonemptyMS/MSIDs andtargetprob>=0.75',secondary_rule='MULTI-SECPEP retainedseparately,notmergedinto primarypass',MBR='MULTI-MATCH retainedseparately;not standaloneMSMSproof',no_new_effect_test=True,original_null_retained=True))
evidence=read(src);sites={x['Amino acid']+x['Position']:x for x in read(D/'SOURCE_TARGET_ROWS.json')};obs=read(D/'OBSERVATION_AUDIT.json')
human=''.join(x.strip() for x in (P/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not x.startswith('>'))
records=[];details=[];issues=[]
for ob in obs:
 site=ob['site'];pos=int(site[1:]);sr=sites[site];sid=sr['id'];eids=set(sr['Evidence IDs'].split(';'));candidates=[]
 for e in evidence:
  if e['Experiment']!=ob['experiment']:continue
  a=e['id'] in eids;b=sid in e['Phospho (STY) site IDs'].split(';')
  if a!=b:issues.append(dict(case=ob['case'],site=site,evidence_id=e['id'],site_to_evidence=a,evidence_to_site=b))
  if not(a and b):continue
  seq=e['Sequence'];starts=[m.start() for m in re.finditer('(?='+re.escape(seq)+')',human)]
  probs=e['Phospho (STY) Probabilities'];tokens=re.findall(r'([A-Z])(?:\(([0-9.eE+\-]+)\))?',probs);parsed=''.join(x[0] for x in tokens)
  prob=None;reason=None
  if len(starts)!=1:reason='nonunique_or_unmapped_sequence'
  elif parsed!=seq:reason='missing_or_unparsed_probability_sequence'
  else:
   offset=pos-starts[0]-1
   if not(0<=offset<len(tokens)):reason='target_outside_sequence'
   else:prob=float(tokens[offset][1]) if tokens[offset][1] else 0.0
  msms=float(e['MS/MS Count'] or '0')>0 and bool(e['MS/MS IDs'])
  direct=e['Type'] in ['MSMS','MULTI-MSMS'] and msms and prob is not None and prob>=.75
  x=dict(case=ob['case'],experiment=ob['experiment'],site=site,evidence_id=e['id'],type=e['Type'],sequence=seq,probability_string=probs,target_probability=prob,has_msms=msms,direct_localization_pass=direct,parse_issue=reason,proteins=e['Proteins'],MSMS_IDs=e['MS/MS IDs'])
  candidates.append(x);details.append(x)
 records.append(dict(case=ob['case'],experiment=ob['experiment'],site=site,linked_evidence=len(candidates),type_counts=dict(collections.Counter(x['type'] for x in candidates)),direct_support_count=sum(x['direct_localization_pass'] for x in candidates),any_direct_localization_pass=any(x['direct_localization_pass'] for x in candidates),max_probability=max([x['target_probability'] for x in candidates if x['target_probability'] is not None],default=None)))
passing=[c for c in sorted(set(x['case'] for x in records)) if all(x['any_direct_localization_pass'] for x in records if x['case']==c)]
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),measurements=len(records),reciprocal_link_disagreements=len(issues),linked_evidence_rows=len(details),probability_parse_issues=sum(x['parse_issue'] is not None for x in details),direct_localization_supported_measurements=sum(x['any_direct_localization_pass'] for x in records),by_site={s:dict(total=sum(x['site']==s for x in records),direct_supported=sum(x['site']==s and x['any_direct_localization_pass'] for x in records)) for s in ['S283','S297']},all_four_supported_patients=passing,all_four_supported_n=len(passing),new_effect_tests=0)
save(W/'MEASUREMENT_LINKS.json',records);save(W/'EVIDENCE_DETAILS.json',details);save(W/'LINK_DISAGREEMENTS.json',issues);save(W/'RESULTS.json',result);print(json.dumps(result))
