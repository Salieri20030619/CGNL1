from pathlib import Path,PurePosixPath
import csv,json,hashlib,re,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024';A=D/'rollup_audit';A.mkdir(exist_ok=True)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
receipts=json.loads((D/'UMich/UMich_peptide_multisite_receipts.json').read_text())
for r in receipts:assert sha(D/'UMich'/PurePosixPath(r['member']).name)==r['sha256']
save(A/'SPEC_BEFORE_COMPARISON.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),purpose='Identify sourcephosphoforms behindexisting C01;no effectselection',comparisons=['S198single row versusS198S202multi row','S202single row versusS202monophospho row','S297single row versusS297monophospho row'],fields='all shared275exportchannelcolumns includingQC;not275independentbiologicalsamples',no_new_P=True))
tables={};heads={}
for level in ['single-site','multi-site','peptide']:
 f=D/'UMich'/f'Report_abundance_groupby={level}_protNorm=MD_gu=2.tsv'
 with f.open() as h:
  r=csv.DictReader(h,delimiter='\t');heads[level]=r.fieldnames;tables[level]=[x for x in r if x['Gene']=='ENSG00000128849.10']
assert heads['single-site']==heads['multi-site']==heads['peptide']
cols=heads['single-site'][5:];pairs=[]
for single,multi in [('S198','_2_2_S198S202'),('S202','_1_1_S202'),('S297','_1_1_S297')]:
 sr=[x for x in tables['single-site'] if x['Index'].endswith('_'+single)];mr=[x for x in tables['multi-site'] if x['Index'].endswith(multi)];assert len(sr)==len(mr)==1
 a,b=sr[0],mr[0];numeric=[];mismatch=[]
 for c in cols:
  am=a[c]=='NA';bm=b[c]=='NA'
  if am!=bm:mismatch.append(c)
  if not am and not bm:numeric.append(abs(float(a[c])-float(b[c])))
 pairs.append(dict(site=single,single_index=a['Index'],multisite_index=b['Index'],single_peptide=a['Peptide'],multisite_peptide=b['Peptide'],channel_columns=len(cols),joint_numeric=len(numeric),missing_mismatch=mismatch,max_value_difference=max(numeric),reference_intensity_difference=float(a['ReferenceIntensity'])-float(b['ReferenceIntensity'])))
hits=[x for x in tables['multi-site'] if re.search(r'S198(?!\d)',x['Index'].rsplit('_',1)[-1])]
assert len(hits)==1 and hits[0]['Index'].endswith('_2_2_S198S202')
save(A/'TARGET_SOURCE_ROWS.json',tables)
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_row_counts={k:len(v) for k,v in tables.items()},comparisons=pairs,only_exported_multisite_row_localizing_S198=hits[0]['Index'],
 consequence='C01 S198 label represents S198/S202 double-phosphorylated form in inspected official export. Do not claim isolatedS198 occupancy or directfunctional specificity.',
 limitations=['Row equivalence does not replace individual spectra/localization probability','No exportedS198-onlyrow is not proofS198-only neverexistsbiologically','Samecohortprocessingnotindependentreplication','Currentmethodpaper clarifieslogic but old2.4tag hasREADMEonly;exactsourcecodeversionnotestablished'],new_tests=0)
save(A/'RESULTS.json',z);print(json.dumps(z))
