from pathlib import Path
import csv,json,re,hashlib,sys,datetime
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');F=O/'C04_full_fraction_qualification';P=O/'C04_search_pilot'
spec=json.loads((F/'JOINT_IDENTIFICATION_SPEC_BEFORE_RESULTS.json').read_text())
db=P/'mouse_plus_contaminants.fasta';assert hashlib.sha256(db.read_bytes()).hexdigest()==spec['database_sha256']
proteins={};header=None;seq=[]
def commit():
 if header is not None:
  token=header.split()[0];g=re.search(r'\bGN=(\S+)',header);proteins[token]={'sequence':''.join(seq),'gene':g.group(1) if g else None}
for line in db.read_text().splitlines():
 if line.startswith('>'):commit();header=line[1:];seq=[]
 else:seq.append(line.strip())
commit()
canonical=next(v['sequence'] for k,v in proteins.items() if k.startswith('sp|Q6AW69|'))
def positions(sequence,peptide):return [m.start()+1 for m in re.finditer('(?='+re.escape(peptide)+')',sequence)]
def annotate(x,has_q):
 peptide=x['plain_peptide'];mods=[]
 for item in x['modifications'].split(','):
  if not item:continue
  fields=item.split('_');assert len(fields)>=3,item
  if abs(float(fields[2])-79.966331)<.00001:
   pos=int(fields[0]);assert 1<=pos<=len(peptide) and peptide[pos-1] in 'STY';mods.append(pos)
 assert len(set(mods))==len(mods)
 scoretext=x.get('ascore_sitescores','');scorematches=re.findall(r'(\d+):([-+0-9.eE]+)',scoretext);scores={int(k):float(v) for k,v in scorematches}
 source=[]
 for token in x['protein'].split(','):
  assert token in proteins,token
  q=proteins[token];starts=positions(q['sequence'],peptide);assert starts,(token,peptide)
  source.append({'protein':token,'gene':q['gene'],'peptide_starts_1based':starts,'engine_proposed_sites':[[{'position':s+i-1,'residue':peptide[i-1]} for i in mods] for s in starts]})
 starts=positions(canonical,peptide);acceptors=[i+1 for i,a in enumerate(peptide) if a in 'STY']
 sole=bool(mods) and len(mods)==1 and len(acceptors)==1
 score_positions_match=set(scores)==set(mods) if mods else not scores
 provisional=bool(mods) and score_positions_match and all(v>=13 for v in scores.values())
 return {'file':x.get('file','pilot_fraction1'),'mode':x.get('mode','HCD'),'scan':int(x['scan']),'charge':int(x['charge']),'peptide':peptide,'phospho_count':len(mods),'engine_proposed_peptide_positions':mods,'engine_proposed_canonical_sites':[[s+i-1 for i in mods] for s in starts],'canonical_mapping_unique':len(starts)==1,'protein_mappings':source,'all_reported_proteins_Cgnl1':all(z['gene'] and z['gene'].lower()=='cgnl1' for z in source),'eligible_acceptors':acceptors,'sole_acceptor_within_peptide':sole,'site_scores':scores,'site_score_positions_match_modifications':score_positions_match,'site_scores_ge13_provisional_only':provisional,'MOB_not_site_probability':x.get('ascorepro'),'q_global':x.get('q_global'),'q_phospho':x.get('q_phospho'),'PSM_dual_q_pass':bool(has_q and mods and float(x['q_global'])<=.01 and float(x.get('q_phospho',1))<=.01),'classification':'unambiguous_acceptor_within_peptide' if sole else 'provisional_localization_unvalidated_FLR' if provisional else 'localization_unresolved','HCD_CID_same_precursor_agreement':'pending_mzML_linkage','PSM_not_biological_replicate':True}
if '--pilot-check' in sys.argv:
 with (P/'run_v2_all_spectra/pilot_hcd.txt').open() as h:
  h.readline();rows=[x for x in csv.DictReader(h,delimiter='\t') if '|Q6AW69' in x['protein'] and 'DECOY_' not in x['protein']]
 checks=[annotate(x,False) for x in rows];first=next(x for x in checks if x['scan']==26694)
 assert first['engine_proposed_peptide_positions']==[6] and first['engine_proposed_canonical_sites']==[[262]]
 assert first['classification']=='localization_unresolved' and first['site_scores']=={6:0.0} and not first['sole_acceptor_within_peptide']
 result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pilot_records_checked':len(checks),'static_TMT_not_misclassified_as_phosphorylation':True,'known_unresolved_site_retained':True,'full_scope_annotation_not_run':True};(F/'ANNOTATOR_READINESS.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result));raise SystemExit(0)
A=F/'joint_identification'
if not (A/'RESULTS.json').exists():raise SystemExit('Joint inputs pending; no candidate annotation produced')
manifest=json.loads((A/'INPUT_MANIFEST.json').read_text());assert manifest['complete_frozen_file_scope']
result=json.loads((A/'RESULTS.json').read_text())
with (A/'INDEPENDENT_R_QA.csv').open() as h:qa=next(csv.DictReader(h))
assert int(qa['rows'])==result['total_PSMs'] and float(qa['global_max_error'])<1e-12 and float(qa['phospho_max_error'])<1e-12
assert (A/'INDEPENDENT_R_QA.csv').stat().st_mtime>=(A/'PSM_QVALUES.tsv.gz').stat().st_mtime
rows=json.loads((A/'CGNL1_ALL_CANDIDATES.json').read_text());annotated=[annotate(x,True) for x in rows]
for row in annotated:
 if row['phospho_count']==0:row['classification']='unphosphorylated_PSM'
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'joint_input_manifest_sha256':hashlib.sha256((A/'INPUT_MANIFEST.json').read_bytes()).hexdigest(),'candidates':annotated,'new_biological_tests':0,'site_FDR_not_established':True,'precursor_linkage_quantification_and_protein_gates_pending':True,'whole_project_complete':False}
(A/'CGNL1_SITE_ANNOTATIONS.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8');print('Annotated',len(annotated),'PSMs; no disease inference')
