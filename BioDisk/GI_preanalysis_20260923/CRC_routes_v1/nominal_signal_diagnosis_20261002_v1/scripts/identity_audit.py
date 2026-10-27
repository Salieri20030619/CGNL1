from pathlib import Path
import json,gzip,datetime
import numpy as np,pandas as pd
from scipy import io
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parents[1];B=R.parent;L=B/'threshold_ladder_20261002_v1';D=R/'results';S=json.loads((R/'IDENTITY_AUDIT_SPECIFICATION.json').read_text());C=json.loads((L/'SPECIFICATION.json').read_text());cache=json.loads((L/'CACHE_INDEX.json').read_text());idx={(x['dataset'],x['sample']):x for x in cache};rows=[];ident=[]
I=B/'spatial_v3/background_validation/run_20260924T052010Z';jobs=[]
for s in json.loads((I/'SAMPLES.json').read_text()):
 if s['dataset']=='GSE285991':jobs.append((s['dataset'],s['sample'],[I/x for x in s['files']]))
for s in json.loads((B/'cohort_verification_20260924/GSE245552_LM_SAMPLES.json').read_text()):jobs.append(('GSE245552',s['sample'],[B/'cohort_verification_20260924/GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]))
for ds,sample,paths in jobs:
 rec=idx[(ds,sample)];fs={k:next(p for p in paths if k in p.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(fs['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(fs['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
 with gzip.open(fs['matrix'],'rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsr()
 with np.load(rec['path']) as f:z={k:f[k] for k in f.files}
 ix=pd.Index(bc).get_indexer(z['barcode']);assert (ix>=0).all();X=X[:,ix];assert np.allclose(np.asarray(X.sum(axis=0)).ravel(),z['lib'])
 def vals(gs):
  out=[]
  for g in gs:
   ii=np.flatnonzero(sy==g);out.append(np.asarray(X[ii].sum(axis=0)).ravel() if len(ii) else np.full(len(ix),np.nan))
  return np.array(out)
 ec=vals(S['support_genes']);original=vals(C['identity_markers']['endothelial']);assert np.isfinite(ec).all();q=(z['qlib']>=500)&(z['qng']>=200)&(z['qmt']<=25)
 for group,k in [('one_marker_removed',q&(z['hit']==1)),('two_plus_retained',q&(z['hit']>=2))]:
  n=int(k.sum())
  if not n:continue
  support=(ec[:,k]>0).sum(axis=0);row=dict(dataset=ds,donor=rec['donor'],sample=sample,group=group,cells=n,median_UMI=float(np.median(z['lib'][k])),fraction_independent_support_ge1=float((support>=1).mean()),fraction_independent_support_ge2=float((support>=2).mean()),mean_independent_endo_score=float(np.log1p(ec[:,k]/z['lib'][k]*1e4).mean()))
  for gene,v in zip(S['support_genes'],ec):row[gene+'_fraction']=float((v[k]>0).mean())
  for lineage in ['immune','epithelial','fibroblast','hepatocyte']:
   v=vals(C['identity_markers'][lineage])[:,k];row[lineage+'_score']=float(np.nanmean(np.log1p(v/z['lib'][k]*1e4)))
  rows.append(row)
  if group=='one_marker_removed':
   sole=(original[:,k]>0).argmax(axis=0)
   for j,gene in enumerate(C['identity_markers']['endothelial']):ident.append(dict(dataset=ds,donor=rec['donor'],sample=sample,sole_marker=gene,cells=int((sole==j).sum()),fraction=float((sole==j).mean())))
 print('IDENTITY_AUDIT',ds,sample,flush=True)
t=pd.DataFrame(rows);t.to_csv(D/'independent_identity_lesion_profiles.tsv',sep='\t',index=False);p=t.groupby(['dataset','donor','group'],as_index=False).mean(numeric_only=True);p.to_csv(D/'independent_identity_patient_profiles.tsv',sep='\t',index=False);p.groupby(['dataset','group'],as_index=False).median(numeric_only=True).to_csv(D/'independent_identity_cohort_medians.tsv',sep='\t',index=False);pd.DataFrame(ident).to_csv(D/'one_marker_identity_distribution.tsv',sep='\t',index=False)
(R/'IDENTITY_COMPLETE.json').write_text(json.dumps({'exit_code':0,'libraries':len(jobs),'descriptive_only':True,'no_new_exclusion':True,'finished':datetime.datetime.now().astimezone().isoformat()},indent=2))
