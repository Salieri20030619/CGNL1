from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib,datetime,gzip,re
import numpy as np,pandas as pd,anndata as ad
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'increment_context_20261012_v1';O=R/'results'
for d in [O,R/'figures',R/'logs',R/'metadata']:d.mkdir(parents=True,exist_ok=True)
prior=B/'epithelial_independence_20261012_v1';act=B/'actomyosin_investigation_20261012_v1'
old=json.loads((prior/'SPECIFICATION.json').read_text())
S=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),target='CGNL1',
question='What conditions preserve the epithelial information increment, and what explains its failure to transfer?',
prior_results_known=True,new_target_results_not_yet_computed=True,
scope='Existing matrices and frozen cell definitions only; no target switch; all branches retained; no clinical or causal claim.',
panels=old['panels'],outcomes=old['outcomes'],ordinary_markers=old['ordinary_marker_comparators'],
primary_outcome='actomyosin',base_covariates=['depth','epithelial_identity','mature_colonocyte'],
modifiers=['mature_colonocyte','regenerative_proxy','cycling','hypoxia'],
measurement_variants=['legacy','pooled_allgenes','pooled_common_genes','within_cohort_rank_legacy'],
measurement='Retain original library QC >=20 epithelial cells. Sensitivities sum eligible libraries per patient; original whole-transcriptome denominator or common-measured-gene denominator. External main definition author epithelial + frozen QC; independent identity intersection separately retained.',
models='M1 controls depth, epithelial identity, maturity, and discovery study. Single-cohort slopes HC3 t; pooled heterogeneity allows cohort-specific nuisance slopes and target slopes. State interaction models allow cohort-specific target slopes and nuisance slopes; common target*within-cohort standardized state interaction.',
therapy='Naive/Treated effect modification in15 external patients. Each subgroup has7/8 patients: exploratory only; descriptive unadjusted slopes; do not fit underpowered subgroup multivariable tests. Combined interaction controls M1 and therapy main effect.',
prediction='Patient leave-one-out and leave-one-cohort-out; compare same nuisance baseline to added CGNL1, and to CGNL1+one predefined interaction. All centering/scaling and fitted coefficients from training only. No tuning or best-state selection. No CV-derived inferential P.',
FDR={'measurement_slopes':'4 scopes within each definition/measurement; all reported',
'heterogeneity':'2 planned slope-heterogeneity tests within each definition/measurement',
'state_interactions':'4 prespecified state modifiers in common-gene main definition',
'therapy':'single prespecified interaction, no correction required',
'benchmark':'27 prespecified predictors within each cohort; same primary actomyosin outcome',
'expression_state_links':'4 state modifiers within each cohort, conditional CGNL1 slope'},
inference='Biological unit=patient. HC3 t; full rank; residual df>=5; no invalid P. Nominal P<0.05 candidate hard rule with family FDR retained; n<10 exploratory. Rank sensitivity is association only, not out-of-sample prediction.',
interpretation='Dataset and platform are confounded; no attribution of heterogeneity to platform without matched material. Interactions are exploratory and require cross-cohort predictive support.',
known_limitations=['external data are noninteger decontX-adjusted counts, no raw layer','discovery treatment metadata may be incomplete','same legacy cohorts already used in target hypothesis generation'])
sp=R/'SPECIFICATION.json'
if not sp.exists():
 sp.write_text(json.dumps(S,ensure_ascii=False,indent=2))
 (R/'SPEC_FREEZE.json').write_text(json.dumps(dict(sha256=hashlib.sha256(sp.read_bytes()).hexdigest(),utc=S['created_utc']),indent=2))
else:
 assert hashlib.sha256(sp.read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
 S=json.loads(sp.read_text())
print('SPEC_FROZEN',flush=True)
sources=[prior/'results/A_patient_profiles.tsv',prior/'results/A_library_profiles.tsv',prior/'cache/A_all_gene_pseudobulk.npz',act/'results/D4_patient_profiles.tsv',B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad',prior/'results/A_CONTROL_FREEZE.json']
finger=[]
for p in sources:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(1048576),b''):h.update(v)
 finger.append(dict(path=str(p),size=p.stat().st_size,sha256=h.hexdigest()))
(R/'SOURCE_FINGERPRINTS.json').write_text(json.dumps(finger,indent=2))
assert finger[4]['sha256']==json.loads((B/'formal_v1_20261002/qualification/multiome_matrix_scale.json').read_text())['sha256']
lib=pd.read_csv(sources[1],sep='\t')
cache=np.load(sources[2],allow_pickle=True)
assert np.array_equal(cache['samples'],lib['sample'])
legacy=pd.read_csv(sources[0],sep='\t')
external_old=pd.read_csv(sources[3],sep='\t')
control=json.loads(sources[5].read_text())['controls']
panels={**S['panels'],**S['outcomes']}
needed=sorted(set(['CGNL1']+S['ordinary_markers']+control+sum(panels.values(),[])))
a=ad.read_h5ad(sources[4]);m=a.obs.copy();X=a.X.tocsr().astype(float);g=a.var_names.astype(str)
assert g.is_unique
assert not list(a.layers) and a.raw is None
qc=(m.nCount_RNA.astype(float)>=500)&(m.nFeature_RNA.astype(float)>=200)&(m['percent.mt'].astype(float)<=25)
ident=pd.read_csv(B/'formal_v1_20261002/02_multiome/cell_identity_QC.tsv.gz',sep='\t').set_index('cell_id').reindex(m.index)
assert np.array_equal(qc.to_numpy(),ident.QC.to_numpy())
author=m.Cell_type.astype(str).eq('Epithelial')
common=sorted(set(cache['genes'].astype(str))&set(g))
pd.Series(common,name='gene').to_csv(R/'metadata/common_genes.tsv',sep='\t',index=False)
assert set(sum(panels.values(),[])+['CGNL1'])<=set(common)
cg=pd.Index(cache['genes'].astype(str));ix=cg.get_indexer(common);jx=g.get_indexer(common)
def prof(counts,genes,den):
 vals=np.log1p(np.asarray(counts)/den*1e4)
 dic=dict(zip(genes,vals))
 row={x:dic.get(x,np.nan) for x in needed}
 for name,gs in panels.items():
  vv=[row[x] for x in gs if np.isfinite(row.get(x,np.nan))]
  row[name]=float(np.mean(vv)) if len(vv)>=.8*len(gs) else np.nan
 return row
rows=[];gene_stats=[];qa={}
for (dataset,patient),v in lib[lib.eligible].groupby(['dataset','patient']):
 inds=v.index.to_numpy();counts=cache['counts'][inds].sum(0);total=float(cache['denominators'][inds].sum())
 meta=dict(dataset=dataset,patient=patient,therapy='Unknown',cells=int(v.cells.sum()),libraries=len(v),depth=float(np.average(v.depth,weights=v.cells)),definition='discovery_frozen')
 for variant,den in [('pooled_allgenes',total),('pooled_common_genes',float(counts[ix].sum()))]:
  rows.append(dict(meta,variant=variant,denominator=den,**prof(counts,cg,den)))
 d=legacy[legacy.patient.eq(patient)].iloc[0].to_dict();d.update(variant='legacy',definition='discovery_frozen',therapy='Unknown',libraries=len(v))
 rows.append(d)
 for gene in ['CGNL1']+S['outcomes']['actomyosin']:
  j=cg.get_loc(gene)
  gene_stats.append(dict(dataset=dataset,patient=patient,gene=gene,source='raw_UMI',counts=float(counts[j]),fraction_positive=float(np.average(cache['detection'][inds,j],weights=v.cells)),CP10k=float(counts[j]/total*1e4)))
for definition,mask in [('author_epithelial_QC',author&qc),('frozen_identity_intersection',author&qc&ident.frozen_identity.eq('epithelial'))]:
 for pat in sorted(m.Patient.astype(str).unique()):
  sel=np.flatnonzero((mask&m.Patient.astype(str).eq(pat)).to_numpy())
  if len(sel)<20:continue
  counts=np.asarray(X[sel].sum(0)).ravel();total=float(counts.sum());pid='Ogden2025:'+pat
  therapy=m.iloc[sel].Therapy.astype(str).unique();assert len(therapy)==1
  meta=dict(dataset='Ogden2025',patient=pid,definition=definition,therapy=therapy[0],cells=len(sel),libraries=m.iloc[sel].Sample.nunique(),depth=float(np.log1p(m.iloc[sel].nCount_RNA.astype(float)).mean()))
  for variant,den in [('legacy',total),('pooled_allgenes',total),('pooled_common_genes',float(counts[jx].sum()))]:
   rows.append(dict(meta,variant=variant,denominator=den,**prof(counts,g,den)))
  if definition=='author_epithelial_QC':
   for gene in ['CGNL1']+S['outcomes']['actomyosin']:
    j=g.get_loc(gene)
    gene_stats.append(dict(dataset='Ogden2025',patient=pid,gene=gene,source='noninteger_adjusted',counts=float(counts[j]),fraction_positive=np.nan,CP10k=float(counts[j]/total*1e4)))
profiles=pd.DataFrame(rows)
for definition in ['author_epithelial_QC','frozen_identity_intersection']:
 z=profiles[profiles.definition.eq(definition)&profiles.variant.eq('legacy')].set_index('patient')
 oldv=external_old[external_old.definition.eq(definition)].set_index('patient').reindex(z.index)
 cols=['CGNL1','actomyosin']+S['base_covariates']+S['modifiers']
 qa[definition+'_max_reconstruction_error']=float(np.nanmax(abs(z[cols].to_numpy()-oldv[cols].to_numpy())))
 assert qa[definition+'_max_reconstruction_error']<1e-10
profiles.to_csv(O/'patient_profiles_all_variants.tsv',sep='\t',index=False)
pd.DataFrame(gene_stats).to_csv(O/'target_measurement.tsv',sep='\t',index=False)
m[['Sample','Patient','Therapy','Tissue']].drop_duplicates().to_csv(R/'metadata/multiome_sample_map.tsv',sep='\t',index=False)
# Preserve source clinical metadata and exact available treatment annotations.
clin=[]
paths=list((B/'cohort_verification_20260924/sc_metadata').glob('*'))+[B/'spatial_v3/background_validation/run_20260924T052010Z/GSE285991.soft.gz']
for p in paths:
 if not p.is_file() or p.stat().st_size>20_000_000:continue
 try:
  text=gzip.open(p,'rt',errors='replace').read() if p.suffix=='.gz' else p.read_text(errors='replace')
 except Exception:continue
 if '245552' not in p.name and '285991' not in p.name:continue
 (R/'metadata'/p.name).write_bytes(p.read_bytes())
 for block in re.split(r'(?=\\^SAMPLE = )',text):
  lines=[l for l in block.splitlines() if l.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics','!Sample_source_name','!Sample_library_strategy','!Sample_extract_protocol','!Sample_data_processing'))]
  if lines:clin.append(dict(source=p.name,record='\n'.join(lines)))
pd.DataFrame(clin).to_csv(R/'metadata/clinical_source_excerpts.tsv',sep='\t',index=False)
qa.update(common_genes=len(common),discovery_patients=legacy.patient.nunique(),external_patients=15,external_identity_patients=13,raw_external_layer_available=False,all_pass=True)
(R/'PREPARE_QA.json').write_text(json.dumps(qa,indent=2))
print('PREPARE_COMPLETE',json.dumps(qa),flush=True)
