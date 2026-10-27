from pathlib import Path
import json,hashlib,re,datetime
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_search_pilot'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n');assert json.loads(f.read_text())==z
mouse=W/'mouse_UP000000589.fasta';assert sha(mouse)=='1caa472043709281c4069a253e29a2e9d3691933e0e88b232c1ac7cbb6ee4f49'
c=W/'contaminants_from_GPM_UniProt.fasta';assert sha(c)==json.loads((W/'contaminants_receipt.json').read_text())['sha256']
ms=mouse.read_text();headers=[x.split()[0] for x in ms.splitlines() if x.startswith('>')];assert len(headers)==len(set(headers))==63328;assert any('|Q6AW69|' in x for x in headers)
cs='\n'.join('>CON__'+x[1:] if x.startswith('>') else x for x in c.read_text().splitlines())+'\n'
db=W/'mouse_plus_contaminants.fasta';db.write_text(ms.rstrip()+'\n'+cs)
params={'database_name':str(db),'decoy_search':'1','num_threads':'4','peptide_mass_tolerance_upper':'20.0','peptide_mass_tolerance_lower':'-20.0','peptide_mass_units':'2','isotope_error':'2','search_enzyme_number':'1','allowed_missed_cleavage':'2','variable_mod01':'15.994915 M 0 3 -1 0 0 0.0','variable_mod02':'79.966331 STY 0 3 -1 0 0 97.976896','max_variable_mods_in_peptide':'3','fragment_bin_tol':'0.02','fragment_bin_offset':'0.0','output_txtfile':'1','output_pepxmlfile':'1','output_percolatorfile':'1','num_output_lines':'1','ms_level':'2','activation_method':'HCD','spectrum_batch_size':'2000','max_duplicate_proteins':'-1','peptide_length_range':'6 50','clear_mz_range':'125.0 132.0','add_Nterm_peptide':'229.162932','add_K_lysine':'229.162932','add_C_cysteine':'57.021464'}
s=(W/'comet.params.new').read_text()
for key,val in params.items():
 s,n=re.subn(r'^'+re.escape(key)+r'\s*=.*$',key+' = '+val,s,flags=re.M);assert n==1,key
s=s.replace('[COMET_ENZYME_INFO]','print_ascorepro_score = 2\n\n[COMET_ENZYME_INFO]');f=W/'comet.params';f.write_text(s)
raw=O/'C04_raw_pilot/phospho_fraction_1.mzML';link=W/'pilot_input.mzML'
if not link.exists():link.symlink_to(raw)
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='frozen_before_any_search_results',software='Comet2026.03rev0 fa08489',database_mouse_records=63328,contaminant_records=113,unresolved_contaminant_accessions=['P04188','P04745'],database_sha256=sha(db),parameters_sha256=sha(f),software_sha256=sha(W/'comet.linux.exe'),input_mzML_sha256='adab96e2d5db2be34eedc1dc6529c1799a65a5b72314693a7b65b44c1b609449',parameters=params,
 purpose='Exploratory singlefractionidentificationfeasibility;notexact2014mouseORFSEQUESTreproduction',
 primary_identification='Concatenatedtarget-decoy competition;one top hit per spectrumquery;E-value ascending with ties grouped;estimatedFDR=(decoys+1)/targets;monotonicq via reversecumulative minimum;report q<=.01PSMs andcomplete distributions',
 peptide_and_localization_boundary='PSM-levelq notprotein/siteFDR;AScorePro experimental,do notconverttoPTMProphetprobability orclaim mechanism',
 target_policy='Retain allCGNL1candidate hits before/afterfilter;no target-only search;onefractionnegative notcohortnegative',
 sample_policy='8biologicalmice inTMTchannels,notMS2/MS3counts;no biologicalPfromsearch',
 expansion_gate='Review globalidentificationquality,score/masserror/decoydistribution,localizationsupport andPSM-MS3linkage beforefullfractiondownload',new_disease_tests=0)
save(W/'SEARCH_SPEC_BEFORE_RESULTS.json',spec);print(json.dumps(dict(database_records=63441,params_sha256=sha(f),database_sha256=sha(db),status='frozen')))
