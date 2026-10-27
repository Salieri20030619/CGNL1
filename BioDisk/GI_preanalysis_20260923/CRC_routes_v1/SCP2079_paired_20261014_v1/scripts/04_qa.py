from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,hashlib,itertools,datetime
from scipy import stats,sparse
R=Path(__file__).resolve().parent.parent;O=R/'results';Q=R.parent/'SCP2079_qualification_20261014_v1'
S=json.loads((R/'SPECIFICATION.json').read_text(encoding='utf-8'))
assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text(encoding='utf-8'))['sha256']
D=pd.read_csv(O/'patient_deltas.tsv',sep='\t');M=pd.read_csv(O/'CGNL1_delta_models.tsv',sep='\t');C=pd.read_csv(O/'paired_changes.tsv',sep='\t');P=pd.read_csv(O/'heldout_predictions.tsv',sep='\t');F=pd.read_csv(O/'predictive_increment.tsv',sep='\t')
def subset(measure,pop):
 g=D[D.measurement.eq(measure)]
 if pop.startswith('author'):g=g[g.definition.eq('author_retained')&g.min_cells.ge(20)]
 else:g=g[g.definition.eq('fixed_qc')&g.min_cells.ge(50 if 'ge50' in pop else 20)]
 if pop.endswith('_liver'):g=g[g.tissueSite.eq('Liver')]
 return g.reset_index(drop=True)
err={'slope':0.,'P_HC3':0.,'CI':0.,'FDR':0.,'paired_exact_P':0.,'prediction':0.,'prediction_metric':0.,'profile_reconstruction':0.}
for _,r in M[M.status.eq('estimated')].iterrows():
 g=subset(r.measurement,r.population);features=S['models'][r.model]+['CGNL1']
 X=np.column_stack([np.ones(len(g)),g[features].to_numpy()]);y=g[r.outcome].to_numpy()
 inv=np.linalg.inv(X.T@X);b=inv@X.T@y;res=y-X@b;h=np.einsum('ij,jk,ik->i',X,inv,X)
 V=inv@(X.T@((res/(1-h))[:,None]**2*X))@inv;se=np.sqrt(V[-1,-1]);df=len(y)-X.shape[1]
 p=2*stats.t.sf(abs(b[-1]/se),df);crit=stats.t.ppf(.975,df);ci=[b[-1]-crit*se,b[-1]+crit*se]
 err['slope']=max(err['slope'],abs(b[-1]-r.slope));err['P_HC3']=max(err['P_HC3'],abs(p-r.P));err['CI']=max(err['CI'],max(abs(np.array(ci)-np.array([r.CI_low,r.CI_high]))))
for table,family_size in [(M,2),(C,3)]:
 for family,g in table.groupby('family'):
  t=g.dropna(subset=['P']);p=np.r_[t.P.to_numpy(),np.ones(family_size-len(t))]
  order=np.argsort(p);vals=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);out=np.empty(len(p));out[order]=vals
  if len(t):err['FDR']=max(err['FDR'],float(np.max(np.abs(out[:len(t)]-t.FDR.to_numpy()))))
for _,r in C[C.status.eq('tested')].iterrows():
 d=subset(r.measurement,r.population)[r.readout].to_numpy();d=d[d!=0]
 if len(d):
  ranks=np.rint(2*stats.rankdata(abs(d))).astype(int);s=np.array(list(itertools.product([0,1],repeat=len(d))),dtype=np.int8)@ranks
  obs=ranks[d>0].sum();p=min(1,2*min(float(np.mean(s<=obs)),float(np.mean(s>=obs))))
 else:p=1.
 err['paired_exact_P']=max(err['paired_exact_P'],abs(p-r.P))
for key,gp in P.groupby(['measurement','population','model','outcome']):
 measure,pop,model,outcome=key;g=subset(measure,pop)
 for _,r in gp.iterrows():
  tr=g[g.donor_id.ne(r.donor_id)];te=g[g.donor_id.eq(r.donor_id)];assert len(te)==1
  for augmented,label in [(False,'baseline'),(True,'augmented')]:
   features=S['models'][model]+(['CGNL1'] if augmented else [])
   X=np.column_stack([np.ones(len(tr)),tr[features]]);T=np.r_[1.,te[features].iloc[0].to_numpy()]
   pred=float(T@np.linalg.lstsq(X,tr[outcome],rcond=None)[0]);err['prediction']=max(err['prediction'],abs(pred-r[label]))
 be=(gp.baseline-gp.observed)**2;ae=(gp.augmented-gp.observed)**2
 f=F[F.measurement.eq(measure)&F.population.eq(pop)&F.model.eq(model)&F.outcome.eq(outcome)].iloc[0]
 assert len(gp)==len(g)==f.n and gp.donor_id.is_unique
 err['prediction_metric']=max(err['prediction_metric'],abs(float(1-ae.mean()/be.mean())-f.relative_MSE_improvement))
# Independently reconstruct every gene/profile with dense selected-gene vectors and metadata denominators.
X=sparse.load_npz(O/'matrix_verified.npz');feat=pd.read_csv(O/'matrix_features.tsv',sep='\t',header=None);bar=pd.read_csv(O/'matrix_barcodes.tsv',sep='\t',header=None)[0]
m=pd.read_csv(Q/'metadata_all_cells.tsv',sep='\t').set_index('NAME').loc[bar].reset_index()
profiles=pd.read_csv(O/'patient_time_profiles.tsv',sep='\t')
genes=json.loads((R/'EXTRACTION_COMPLETE.json').read_text(encoding='utf-8'))['genes']
idx=[int(feat.index[feat[1].eq(g)][0]) for g in genes];A=X[idx,:].toarray();del X
qc=(m.total_counts>=500)&(m.nFeature_RNA>=200)&(m.pct_counts_mt<=25)&(m.doublet_score<=.3)
for _,r in profiles.iterrows():
 mask=m.donor_id.eq(r.donor_id)&m.treatment.eq(r.treatment)
 if r.definition=='fixed_qc':mask&=qc
 assert int(mask.sum())==r.cells
 a=A[:,mask.to_numpy()]
 if r.measurement=='pooled_raw_CP10k':v=np.log1p(a.sum(axis=1)/m.loc[mask,'total_counts'].sum()*10000)
 elif r.measurement=='mean_author_log':v=a.mean(axis=1)
 else:v=np.log1p(np.expm1(a).mean(axis=1))
 err['profile_reconstruction']=max(err['profile_reconstruction'],float(np.max(np.abs(v-r[genes].to_numpy(dtype=float)))))
 for key,gs in {**S['controls'],**S['outcomes']}.items():err['profile_reconstruction']=max(err['profile_reconstruction'],abs(float(r[gs].mean())-r[key]))
assert all(v<1e-6 for v in err.values()),err
assert M.loc[M.n.lt(10),'P'].isna().all() and C.loc[C.n.lt(10),'P'].isna().all()
qa={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'max_absolute_errors':err,'spec_unchanged':True,'all_prediction_folds_complete':bool(F.all_folds_complete.all()),'sub10_patient_groups_have_no_inferential_P':True,'independent_profile_reconstruction':True,'numerical_tolerance':1e-6}
(R/'NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa,indent=2),flush=True)