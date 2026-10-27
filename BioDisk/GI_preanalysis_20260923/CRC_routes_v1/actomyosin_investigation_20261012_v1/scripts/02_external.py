from common import *
import anndata as ad
B=R.parent;prior=B/'epithelial_independence_20261012_v1'
panels=json.loads((prior/'SPECIFICATION.json').read_text())['panels'];needed=sorted({'CGNL1'}|set(S['genes'])|set(sum(panels.values(),[])))
source=B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad'
h=hashlib.sha256()
with source.open('rb') as f:
 for block in iter(lambda:f.read(1048576),b''):h.update(block)
oldscale=json.loads((B/'formal_v1_20261002/qualification/multiome_matrix_scale.json').read_text());assert h.hexdigest()==oldscale['sha256']
a=ad.read_h5ad(source);m=a.obs.copy();X=a.X.tocsr().astype(np.float64);g=a.var_names.astype(str);assert g.is_unique
den=np.asarray(X.sum(1)).ravel();assert (den>0).all()
qc=(m.nCount_RNA.astype(float)>=500)&(m.nFeature_RNA.astype(float)>=200)&(m['percent.mt'].astype(float)<=25)
old=pd.read_csv(B/'formal_v1_20261002/02_multiome/cell_identity_QC.tsv.gz',sep='\t').set_index('cell_id').reindex(m.index);assert old.index.equals(m.index) and np.array_equal(old.QC.to_numpy(),qc.to_numpy())
author=m.Cell_type.astype(str).eq('Epithelial');conf=author&old.frozen_identity.eq('epithelial')
m['depth']=np.log1p(m.nCount_RNA.astype(float));m['patient']=m.Patient.astype(str);m['subtype']=m.Cell_subtype.astype(str)
present=[x for x in needed if x in g];idx=g.get_indexer(present);V=X[:,idx]
assert all(x in present for x in S['genes']+['CGNL1'])
def profile(ix):
 total=np.asarray(V[ix].sum(0)).ravel();lib=float(den[ix].sum());val=np.log1p(total/max(lib,1)*1e4);row=dict(zip(present,val))
 for name,gs in panels.items():
  vv=[row[x] for x in gs if x in row];row[name]=float(np.mean(vv)) if len(vv)>=.8*len(gs) else np.nan
 row['actomyosin']=float(np.mean([row[x] for x in S['genes']]))
 row['depth']=float(m.iloc[ix].depth.mean());row['adjusted_library_sum']=lib
 return row
rows=[];bins=[];coverage=[]
for patient in sorted(m.patient.unique()):
 patientmask=m.patient.eq(patient)
 for definition,identity in [('author_epithelial_QC',author),('frozen_identity_intersection',conf)]:
  mask=(patientmask&qc&identity).to_numpy();ix=np.flatnonzero(mask);coverage.append(dict(patient=patient,definition=definition,cells=len(ix),eligible=len(ix)>=20))
  if len(ix)>=20:
   therapies=sorted(m.iloc[ix].Therapy.astype(str).unique());assert len(therapies)==1
   rows.append(dict(patient='Ogden2025:'+patient,definition=definition,dataset='Ogden2025',therapy=therapies[0],cells=len(ix),**profile(ix)))
 if int((patientmask&qc&author).sum())>=20:
  for subtype in sorted(m.loc[patientmask&qc&author,'subtype'].unique()):
   ix=np.flatnonzero((patientmask&qc&author&m.subtype.eq(subtype)).to_numpy())
   if len(ix)>=20:bins.append(dict(patient='Ogden2025:'+patient,subtype=subtype,cells=len(ix),**profile(ix)))
pd.DataFrame(coverage).to_csv(O/'D4_cell_coverage.tsv',sep='\t',index=False)
v=pd.DataFrame(rows);v.to_csv(O/'D4_patient_profiles.tsv',sep='\t',index=False);bins=pd.DataFrame(bins);bins.to_csv(O/'D4_patient_subtype_profiles.tsv',sep='\t',index=False)
models=[];elig=[];prs=[];perfs=[]
for definition,z in v.groupby('definition'):
 for model,cols in S['models'].items():
  fit=fitted(z,cols,'actomyosin')
  if fit is None:elig.append(dict(definition=definition,model=model,status='not_estimable'));continue
  row,u,D=fit;models.append(dict(definition=definition,model=model,family='D4_three_models_'+definition,**row))
  elig.append(dict(definition=definition,model=model,status='estimated',n=row['n']))
  p=crossval(z,cols,'actomyosin');p['definition']=definition;p['model']=model;prs.append(p);perfs.append(dict(definition=definition,model=model,**performance(p)))
 cols=S['models']['M1']+['therapy'];fit=fitted(z,cols,'actomyosin')
 if fit:models.append(dict(definition=definition,model='M1_therapy',family='D4_therapy_single_'+definition,**fit[0]))
t=pd.DataFrame(models);t['FDR']=np.nan
for family,ix in t.groupby('family').groups.items():t.loc[ix,'FDR']=bh(t.loc[ix,'P'],1 if 'single' in family else 3)
t.to_csv(O/'D4_external_models.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(O/'D4_model_eligibility.tsv',sep='\t',index=False)
pd.concat(prs).to_csv(O/'D4_heldout_predictions.tsv',sep='\t',index=False);pd.DataFrame(perfs).to_csv(O/'D4_predictive_increment.tsv',sep='\t',index=False)
# Within-patient and within-author-subtype secondary model.
bb=bins.dropna(subset=['CGNL1','depth','epithelial_identity','mature_colonocyte','actomyosin']).reset_index(drop=True)
G=bb.patient.nunique();basecols=['depth','epithelial_identity','mature_colonocyte']
D=design(bb,basecols)
for col in ['patient','subtype']:
 for level in sorted(bb[col].unique())[1:]:D[col+'_'+level]=bb[col].eq(level).astype(float)
D=D[[x for x in D if x!='CGNL1']+['CGNL1']]
state={'patients':int(G),'bins':len(bb),'inference_unit':'patient_cluster'}
if G>=10 and len(bb)>=30 and np.linalg.matrix_rank(D)==D.shape[1] and len(bb)-D.shape[1]>=5:
 mm=sm.OLS(bb.actomyosin,D).fit(cov_type='cluster',cov_kwds={'groups':bb.patient,'use_correction':True,'df_correction':True},use_t=True)
 state.update(status='estimated',slope=float(mm.params.CGNL1),beta=float(mm.params.CGNL1*bb.CGNL1.std()/bb.actomyosin.std()),P=float(mm.pvalues.CGNL1),FDR=float(mm.pvalues.CGNL1),family='D4_within_subtype_single_target',effective_inference_df=int(G-1))
else:state['status']='rank_or_coverage_failure'
(O/'D4_within_subtype.json').write_text(json.dumps(state,indent=2))
# Cross-platform transfer is reported without any target-cohort fitting.
discovery=pd.read_csv(prior/'results/A_model_input.tsv',sep='\t');external=v[v.definition.eq('author_epithelial_QC')]
cols=[x for x in S['models']['M1'] if x!='study'];trs=[]
for label,tr,te in [('discovery_to_multiome',discovery,external),('multiome_to_discovery',external,discovery)]:
 bpred=predict_raw(tr,te,cols,'actomyosin',False);apred=predict_raw(tr,te,cols,'actomyosin',True)
 if bpred is None or apred is None:continue
 for i,(_,rr) in enumerate(te.iterrows()):trs.append(dict(direction=label,patient=rr.patient,observed=rr.actomyosin,baseline=float(bpred[i]),augmented=float(apred[i])))
pd.DataFrame(trs).to_csv(O/'D4_transfer_predictions.tsv',sep='\t',index=False)
(O/'D4_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':h.hexdigest(),'source_scale':'author_decontX_adjusted_counts','models':len(t),'source_patients':int(m.patient.nunique()),'author_epithelial_QC_cells':int((author&qc).sum()),'frozen_intersection_cells':int((conf&qc).sum()),'raw_detection_P_not_used':True},indent=2))
print('D4_COMPLETE',len(t),flush=True)
