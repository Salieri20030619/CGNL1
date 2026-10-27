from pathlib import Path
import json,hashlib,datetime,csv,collections,html
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parent
M=pd.read_csv(R/'metadata_all_cells.tsv',sep='\t',dtype=str,keep_default_na=False)
for col in ['total_counts','nCount_RNA','nFeature_RNA','pct_counts_mt','doublet_score','PFSmo']:M[col]=pd.to_numeric(M[col],errors='raise')
assert M.NAME.is_unique and len(M)==95421
invariant={c:int(M.groupby('donor_id')[c].nunique().max()) for c in ['patientID','tissueSite','MMRstatus','sex','PFSmo','response','RECIST']}
assert all(v==1 for v in invariant.values()),invariant
qc=(M.total_counts>=500)&(M.nFeature_RNA>=200)&(M.pct_counts_mt<=25)&(M.doublet_score<=.3)
M['fixed_qc']=qc
site=M.groupby('donor_id').tissueSite.first()
rows=[]
for label,mask in [('author_retained',np.ones(len(M),dtype=bool)),('fixed_qc',qc)]:
 sub=M.loc[mask]
 c=sub.groupby(['donor_id','treatment']).size().unstack(fill_value=0).reindex(index=site.index,columns=['Pre','On'],fill_value=0).fillna(0).astype(int)
 for donor,row in c.iterrows():
  rows.append({'definition':label,'donor_id':donor,'site':site[donor],'Pre_cells':row['Pre'],'On_cells':row['On'],'paired_any':bool(row.min()>0),'paired_ge20':bool(row.min()>=20),'paired_ge50':bool(row.min()>=50)})
coverage=pd.DataFrame(rows)
coverage.to_csv(R/'patient_pair_eligibility.tsv',sep='\t',index=False,encoding='utf-8')
M.groupby(['donor_id','treatment','tissueSite','sampleID']).agg(n_cells=('NAME','size'),fixed_qc_cells=('fixed_qc','sum')).to_csv(R/'technical_library_coverage.tsv',sep='\t',encoding='utf-8')
clinical=M.groupby('donor_id')[['tissueSite','MMRstatus','sex','PFSmo','response','RECIST']].first()
clinical.to_csv(R/'patient_metadata.tsv',sep='\t',encoding='utf-8')
summary=[]
for label,g in coverage.groupby('definition',sort=False):
 for site_label,h in [('all_sites',g)]+list(g.groupby('site',sort=True)):
  summary.append({'definition':label,'site':site_label,'donors':len(h),'paired_any':int(h.paired_any.sum()),'paired_ge20':int(h.paired_ge20.sum()),'paired_ge50':int(h.paired_ge50.sum())})
pd.DataFrame(summary).to_csv(R/'qualification_summary.tsv',sep='\t',index=False,encoding='utf-8')
print(pd.DataFrame(summary).to_string(index=False),flush=True)
# Independent coverage reconstruction from source exports, without the merged table.
source={}
for col in ['donor_id','treatment','tissueSite']:
 with (R/'metadata'/f'{col}.tsv').open(encoding='utf-8',newline='') as f:source[col]={x['NAME']:x[col] for x in csv.DictReader(f,delimiter='\t')}
ind=collections.Counter((source['donor_id'][k],source['treatment'][k],source['tissueSite'][k]) for k in source['donor_id'])
assert set(source['donor_id'])==set(M.NAME)
assert ind==collections.Counter({k:int(v) for k,v in M.groupby(['donor_id','treatment','tissueSite']).size().items()})
fetch=json.loads((R/'METADATA_FETCH.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((R/'metadata'/(d['name']+'.tsv')).read_bytes()).hexdigest()==d['sha256'] for d in fetch)
qa={'cells':len(M),'unique_cells':M.NAME.nunique(),'metadata_columns':len(fetch),'identical_cell_sets':True,'independent_coverage_reconstruction_match':True,'source_hashes_match':True,'patient_fields_invariant':invariant,'donor_patient_match':bool((M.donor_id==M.patientID).all()),'specimen_biosample_match':bool((M.specimenID==M.biosample_id).all()),'qc_retained_cells':int(qc.sum()),'qc_excluded_cells':int((~qc).sum()),'total_counts_integer_fraction':float(np.isclose(M.total_counts,np.round(M.total_counts),rtol=0,atol=1e-8).mean()),'nCount_RNA_integer_fraction':float(np.isclose(M.nCount_RNA,np.round(M.nCount_RNA),rtol=0,atol=1e-8).mean()),'all_metadata_checks_pass':True,'matrix_checks_complete':False}
(R/'METADATA_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
prior=json.loads(Path(r'C:\Users\31937\CGNL1_new_cohorts_20261013_v1\delivery\SPECIFICATION.json').read_text(encoding='utf-8'))
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'preanalysis_conditional_plan; final measurement freeze awaits matrix qualification','dataset':'SCP2079','target_expression_or_association_results_seen':False,'biological_unit':'donor_id; merge technical libraries within donor/time point; never use cells or library labels as independent n','primary_population':'all author epithelial cells passing fixed QC; donor with >=20 cells at both Pre and On','sensitivity_population':'>=50 cells at both times; retain author-QC-only results as sensitivity','liver_role':'only tissueSite=Liver; if n<10 descriptive only per preceding validation protocol; do not overwrite Other/Muscle/Lymph Node with ontology liver','qc':{'total_counts_min':500,'nFeature_RNA_min':200,'pct_counts_mt_max':25,'doublet_score_max':.3,'note':'metadata QC provisional until matrix/read-count reconciliation; nCount_RNA is fractional and not UMI depth'},'target':'CGNL1','outcomes':prior['outcomes'],'controls':{'epithelial_identity':prior['panels']['epithelial_identity'],'mature_colonocyte':prior['panels']['mature_colonocyte']},'minimum_gene_coverage':'all 7 actomyosin genes; >=80% of each other fixed program; target present','matrix_gate':'Check dimensions, gzip integrity, SHA256, unique features/barcodes, full metadata alignment, nonnegative finite values, integer status, column sums vs total_counts and nCount_RNA. counts filename alone does not prove raw UMI. Do not reverse-normalize into invented raw counts.','measurement_if_verified_raw_counts':'Sum counts within patient/time, ln(1+CP10k) over all measured genes; arithmetic mean of fixed genes per program; depth=log mean cell UMI; delta=On-Pre','measurement_if_only_normalized':'Stop and document units; freeze explicit author-scale patient mean strategy separately before viewing target outcomes; no raw-count pseudobulk or inherited raw-count coefficients','primary_models':'Delta outcome ~ delta CGNL1 + delta depth + delta epithelial_identity + delta mature_colonocyte. HC3 t; full rank; residual df>=5; >=10 independent paired donors. M0 depth-only sensitivity. Preserve all estimates/P/FDR.','multiplicity':'BH across the two outcome hypotheses within each population/model/measurement family; nominal gate target P<0.05; label FDR results separately','secondary_paired_changes':'CGNL1, actomyosin, cycling On-Pre, signed effect and patient consistency; separate three-readout family; only infer if valid independent repeats','prediction':'LOPO delta outcome baseline controls vs controls+CGNL1; all scaling fitted in training patients; require all folds estimable; report relative MSE and improved patient counts; no tuning','interpretation':'Single-arm BRAF/MEK/PD-1 treatment context, not CGNL1 perturbation, not isolated drug causal effect, not a general liver-colonization validation','clinical_outcomes':'PFS/RECIST descriptive metadata only in this round; no opportunistic survival screening'}
(R/'CONDITIONAL_ANALYSIS_PLAN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
portal=json.loads((R/'sources/site_api.txt').read_text(encoding='utf-8'))
required=[x for x in portal['study_files'] if x['bucket_location'].startswith('expression_epi_counts/')]
(R/'REQUIRED_DOWNLOADS.json').write_text(json.dumps({'authentication_required':True,'files':required,'optional_normalized_reference':[x for x in portal['study_files'] if x['bucket_location'].startswith('expression_epi_norm/')]},indent=2),encoding='utf-8')
status={'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':'SCP2079 qualification before CGNL1 paired analysis','metadata_qualification_complete':True,'matrix_scale_qualification_complete':False,'inferential_analysis_started':False,'whole_project_complete':False,'download_block':'Official download HTTP401; anonymous published GCS object HTTP403; authenticated user download required','summary':summary,'next_action':'Obtain expression_epi_counts matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz via authorized portal sign-in; qualify actual scale before final analysis freeze','server_directory':'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/SCP2079_qualification_20261014_v1'}
(R/'PROGRESS_SNAPSHOT.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2),flush=True)