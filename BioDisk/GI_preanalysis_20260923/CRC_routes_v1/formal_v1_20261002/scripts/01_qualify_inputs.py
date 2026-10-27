from pathlib import Path
import json,gzip,hashlib,datetime,sys,platform,collections,importlib.metadata
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=Path(__file__).resolve().parents[1];C=json.loads((R/'frozen/EXECUTION_CONFIG.json').read_text())
Q=R/'qualification';D=R/'registry'
def save(n,x):
 (R/n).write_text(json.dumps(x,indent=2,default=str))
def digest(f):
 h=hashlib.sha256()
 with open(f,'rb') as s:
  for z in iter(lambda:s.read(8*1024**2),b''):h.update(z)
 return h.hexdigest()
def tab(df,p):df.to_csv(p,sep='\t',index=False)
# Reuse prior integrity records and preserve their provenance, do not rehash 202 GB.
inputs=[]
for wave in ['wave01_core','wave02_extended','wave03_supplements']:
 mf=B/'data_prepare_20260929'/wave/'MANIFEST.json'
 for x in json.loads(mf.read_text()):
  f=Path(x['path']);assert f.is_file() and f.stat().st_size==x['bytes']
  inputs.append(dict(dataset=x['dataset'],path=str(f),bytes=x['bytes'],sha256=x['sha256'],hash_source=str(mf),integrity='prior_verified_current_size_checked'))
p=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_data_prepare_20260929')
assert json.loads((p/'COMPLETE.json').read_text())['verified']==304
for x in json.loads((p/'MANIFEST.json').read_text()):
 f=Path(x['path']);assert f.is_file() and f.stat().st_size==x['bytes'] and x['status']=='downloaded_verified'
 inputs.append(dict(dataset='E-MTAB-13652_ATAC',path=str(f),bytes=x['bytes'],sha256=x['sha256'],md5=x['md5'],hash_source=str(p/'MANIFEST.json'),integrity='prior_official_md5_verified_current_size_checked'))
tab(pd.DataFrame(inputs),D/'downloaded_input_registry.tsv')
old=json.loads((B/'data_prepare_20260929/EXISTING_DATA_CATALOG.json').read_text())
for x in old:x['exists_now']=Path(x['path']).exists()
save('registry/reused_input_catalog.json',old)
# Parse official GEO sample information; title is a sample label, not a confirmed patient.
meta={};gsm=None
with gzip.open(B/'expanded_screen_20260925/metadata/GSE234804.soft.gz','rt') as f:
 for l in f:
  if l.startswith('^SAMPLE = '):gsm=l.strip().split(' = ',1)[1];meta[gsm]={}
  elif gsm and l.startswith(('!Sample_title = ','!Sample_source_name_ch1 = ')):
   k,v=l.strip().split(' = ',1);meta[gsm][k]=v
rows=[];lineage_rows=[];coverage=[];sc_registry=[]
labels=list(C['lineage_markers']);panel=C['junction_panel'];target=C['target']
for f in sorted((B/'data_prepare_20260929/wave01_core/GSE234804').glob('*.h5seurat')):
 gsm=f.name.split('_')[0];title=meta[gsm]['!Sample_title'];tissue=meta[gsm]['!Sample_source_name_ch1']
 sc_registry.append(dict(dataset='GSE234804',sample_id=gsm,sample_label=title,tissue=tissue,patient_id='',patient_status='unresolved_do_not_count_sample_as_patient',independence='pending_original_study_and_cross_cohort_review',path=str(f)))
 with h5py.File(f) as h:
  genes=h['assays/RNA/features'].asstr()[:];cells=h['cell.names'].asstr()[:];z=h['assays/RNA/counts'];assert z.shape==(len(cells),len(genes));assert len(set(genes))==len(genes) and len(set(cells))==len(cells)
  pieces=[]
  for start in range(0,len(cells),256):pieces.append(sparse.csr_matrix(z[start:start+256]))
  X=sparse.vstack(pieces,format='csr');X.eliminate_zeros()
 assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data),atol=1e-6,rtol=0)
 lib=np.asarray(X.sum(axis=1)).ravel();ng=X.getnnz(axis=1);mt=np.asarray(X[:,np.char.startswith(genes.astype(str),'MT-')].sum(axis=1)).ravel()/np.maximum(lib,1)*100
 qc=(lib>=500)&(ng>=200)&(mt<=25);scores=[];hits=[]
 for label in labels:
  gs=C['lineage_markers'][label];ix=np.flatnonzero(np.isin(genes,gs));assert len(ix)>0
  small=X[:,ix].toarray();scores.append(np.log1p(small/np.maximum(lib[:,None],1)*1e4).mean(axis=1));hits.append((small>0).sum(axis=1))
 scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2]
 valid=qc&(margin>=.4)&(hits[np.arange(len(cells)),best]>=2)
 assigned=np.where(valid,np.array(labels)[best],'unresolved');assigned[~qc]='QC_failed'
 endo=np.flatnonzero(assigned=='endothelial');cgix=np.flatnonzero(genes==target);assert len(cgix)==1;cg=X[:,cgix[0]].toarray().ravel()
 feat=np.c_[np.log1p(lib),np.log1p(ng)];pos=endo[cg[endo]>0];neg=endo[cg[endo]==0];used=set();pairs=[]
 for i in sorted(pos,key=lambda k:(cg[k],k)):
  cand=np.array([k for k in neg if k not in used],dtype=int)
  if not len(cand):break
  delta=np.abs(feat[cand]-feat[i]);ok=(delta<=.25).all(axis=1)
  if ok.any():
   cand=cand[ok];delta=delta[ok];j=cand[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((i,j))
 is_lm='liver metastasis' in tissue.lower()
 rows.append(dict(sample_id=gsm,sample_label=title,tissue=tissue,is_liver_metastasis=is_lm,cells=len(cells),qc_cells=int(qc.sum()),endo_cells=len(endo),CGNL1_positive_endo=len(pos),matched_pairs=len(pairs),passes_cell_coverage50=len(endo)>=50,passes_match_coverage10=len(endo)>=50 and len(pairs)>=10,patient_independence='pending',eligible_confirmatory=False))
 for label in labels+['unresolved','QC_failed']:lineage_rows.append(dict(sample_id=gsm,lineage=label,cells=int((assigned==label).sum())))
 for gene in [target]+panel:coverage.append(dict(dataset='GSE234804',sample_id=gsm,gene=gene,present=bool(np.any(genes==gene))))
 out=Q/gsm;out.mkdir();sparse.save_npz(out/'endo_counts_cell_by_gene.npz',X[endo]);tab(pd.DataFrame({'gene':genes}),out/'genes.tsv')
 tab(pd.DataFrame({'cell':cells,'UMI':lib,'genes':ng,'mito_pct':mt,'QC':qc,'frozen_lineage':assigned,'margin':margin}),out/'cell_qc_lineage.tsv.gz')
 tab(pd.DataFrame({'cell':cells[endo],'UMI':lib[endo],'genes':ng[endo],'CGNL1_count':cg[endo]}),out/'endo_cells.tsv')
 tab(pd.DataFrame(pairs,columns=['positive_row','negative_row']),out/'matched_cell_rows.tsv')
 tab(pd.DataFrame(rows),Q/'GSE234804_library_eligibility.tsv');print('QUALIFIED_LIBRARY',gsm,len(endo),len(pairs),flush=True)
tab(pd.DataFrame(lineage_rows),Q/'GSE234804_lineage_counts.tsv');tab(pd.DataFrame(coverage),Q/'panel_coverage.tsv')
# Fractional decontaminated matrix is qualified structurally but held for raw-detection analyses.
f=B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad'
a=ad.read_h5ad(f);assert a.obs_names.is_unique and a.var_names.is_unique
x=a.X.data if sparse.issparse(a.X) else np.asarray(a.X).ravel();assert np.isfinite(x).all() and (x>=0).all()
scale={'file':str(f),'shape':a.shape,'sha256':digest(f),'noninteger_values':int((abs(x-np.round(x))>1e-6).sum()),'positive_values_below_1e_minus_6':int(((x>0)&(x<1e-6)).sum()),'minimum':float(x.min()),'maximum':float(x.max()),'matrix_role':'decontaminated_expression_not_verified_raw_UMI','raw_detection_analyses':'HOLD until source scale and appropriate raw counts verified','scope':'author cell labels only for initial coverage; independent identity QC remains pending'}
save('qualification/multiome_matrix_scale.json',scale)
obs=a.obs.copy();obs['cell_id']=obs.index
counts=obs.groupby(['Patient','Sample','Therapy','Tissue','Cell_type'],observed=True).size().rename('cells').reset_index();tab(counts,Q/'multiome_author_cell_coverage.tsv')
endo=obs.loc[obs.Cell_type.astype(str)=='Endothelial'];tab(endo,Q/'multiome_author_endothelial_metadata.tsv')
for patient,g in obs.groupby('Patient',observed=True):
 ne=int((g.Cell_type.astype(str)=='Endothelial').sum());sc_registry.append(dict(dataset='CRCLM_multiome',sample_id=';'.join(g.Sample.astype(str).unique()),sample_label=str(patient),patient_id='CRCLM_multiome:'+str(patient),patient_status='author_Patient_field',tissue=';'.join(g.Tissue.astype(str).unique()),independence='RNA_ATAC_same_patient_cross_study_review_pending',path=str(f)))
 print('MULTIOME_AUTHOR_ENDO',patient,ne,flush=True)
tab(pd.DataFrame(sc_registry),D/'new_singlecell_sample_registry.tsv')
# Direct mapping to official ENA run, sample and donor; read/lane counts do not increase n.
sdrf=pd.read_csv(B/'expanded_screen_20260925/metadata/E-MTAB-13652.sdrf.txt',sep='\t',dtype=str)
cols=['Source Name','Characteristics[individual]','Comment[ENA_SAMPLE]','Comment[ENA_RUN]','Comment[SUBMITTED_FILE_NAME]','Comment[READ_TYPE]','Comment[READ_INDEX]']
sdrf=sdrf[cols].drop_duplicates();download=pd.DataFrame(json.loads((p/'MANIFEST.json').read_text()));download['filename']=download['path'].map(lambda v:Path(v).name)
joined=download[['run','filename','path','bytes','md5','sha256']].merge(sdrf,left_on=['run','filename'],right_on=['Comment[ENA_RUN]','Comment[SUBMITTED_FILE_NAME]'],how='left',validate='one_to_one')
assert len(joined)==304;assert joined['Characteristics[individual]'].notna().all()
joined['RNA_patient_present']=joined['Characteristics[individual]'].isin(obs.Patient.astype(str));assert joined.RNA_patient_present.all()
tab(joined,D/'ATAC_RNA_sample_run_mapping.tsv')
# Metadata registries reuse prior reviewed input and retain unresolved fields.
for ds in ['GSE159216','GSE73255']:
 f=B/'cohort_verification_20260924/bulk'/f'{ds}_LM_only_metadata_CGNL1.tsv';m=pd.read_csv(f,sep='\t',dtype=str,keep_default_na=False)
 m=m[[c for c in m if not c.startswith('CGNL1__')]];tab(m,D/f'{ds}_clinical_sample_registry.tsv')
sp=[]
for ds in ['GSE206552','GSE254364','GSE217414']:
 for v in json.loads((B/'expanded_screen_20260925/metadata'/f'{ds}_samples.json').read_text()):sp.append(dict(dataset=ds,sample_id=v['gsm'],title=v['title'],characteristics=json.dumps(v.get('characteristics',[])),patient_id='',overlap_status='pending_NicheDE_original_source_review',geometry_status='pending_physical_scale_and_region_review'))
tab(pd.DataFrame(sp),D/'new_spatial_sample_registry.tsv')
shutil_source=B/'cohort_verification_20260924/qualified_endothelial_patients.tsv'
pd.read_csv(shutil_source,sep='\t').assign(role='previously_explored_not_new_validation',source=str(shutil_source)).to_csv(D/'previously_explored_endo_registry.tsv',sep='\t',index=False)
wet=[{'question':'CGNL1 source and localization','required_evidence':'patient lineage reproducibility plus spatial localization','candidate_status':'pending A03-A12','future_readout':'cell-type-resolved RNA/protein localization','causal_claim':'not established'}, {'question':'endothelial junction state','required_evidence':'A05 association plus A06 depth/subtype controls','candidate_status':'pending','future_readout':'junction organization and barrier function after perturbation','causal_claim':'not established'}, {'question':'upstream regulation','required_evidence':'A17 TF evidence and qualified A18 chromatin support','candidate_status':'pending','future_readout':'candidate regulator perturbation and rescue','causal_claim':'not established'}]
tab(pd.DataFrame(wet),R/'wetlab/validation_bridge.tsv')
summary={'status':'qualification_batch_complete_not_phase1_complete','GSE234804_libraries':len(rows),'GSE234804_LM_libraries':sum(v['is_liver_metastasis'] for v in rows),'GSE234804_LM_libraries_endo50':sum(v['is_liver_metastasis'] and v['passes_cell_coverage50'] for v in rows),'GSE234804_LM_libraries_match10':sum(v['is_liver_metastasis'] and v['passes_match_coverage10'] for v in rows),'new_confirmed_independent_patients':None,'multiome_author_patients':int(obs.Patient.nunique()),'multiome_author_endo_cells':len(endo),'multiome_author_patients_endo50':int((endo.groupby('Patient',observed=True).size()>=50).sum()),'ATAC_files_mapped':len(joined),'ATAC_RNA_patients':int(joined['Characteristics[individual]'].nunique()),'target_association_tests_run':False,'pending':['GSE234804 original-study patient identity and independence','multiome fractional-matrix source and raw counts qualification','A09 spatial original-source overlap and physical geometry','A02 reuse evidence table and final endpoints freeze'],'finished':datetime.datetime.now().astimezone().isoformat()}
save('qualification/QUALIFICATION_SUMMARY.json',summary)
save('ENVIRONMENT.json',{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':{k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','h5py','anndata']}})
outputs=[{'path':str(f.relative_to(R)),'bytes':f.stat().st_size,'sha256':digest(f)} for d in [D,Q,R/'wetlab'] for f in d.rglob('*') if f.is_file()];save('OUTPUT_MANIFEST.json',outputs)
save('BATCH_COMPLETE.json',{'exit_code':0,'scope':'first input qualification batch only','output_files':len(outputs),'finished':summary['finished']});print(json.dumps(summary),flush=True)
