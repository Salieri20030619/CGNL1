from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import json,hashlib,itertools,datetime
import numpy as np,pandas as pd
from scipy.linalg import lstsq
O=Path.cwd();S=json.loads((O/'ANALYSIS_SPECIFICATION.json').read_text());freeze=json.loads((O/'SPEC_FREEZE.json').read_text());d=pd.read_csv(O/'results/spatial_spot_features.tsv.gz',sep='\t');sec=pd.read_csv(O/'results/section_effects.tsv',sep='\t');pt=pd.read_csv(O/'results/patient_effects.tsv',sep='\t');summary=pd.read_csv(O/'results/all_endpoint_summary.tsv',sep='\t');sp=pd.read_csv(O/'results/target_specificity_summary.tsv',sep='\t');covs=S['main_covariates'];targets=S['targets'];errors={};rows=[]
def err(k,a,b):errors[k]=max(errors.get(k,0.),float(np.max(np.abs(np.array(a)-np.array(b)))))
assert hashlib.sha256((O/'ANALYSIS_SPECIFICATION.json').read_bytes()).hexdigest()==freeze['sha256']
# Independently reconstruct eligibility/design and solve with scipy LAPACK QR (gelsy).
for sample,v0 in d.groupby('sample'):
 ref=sec[(sec.variant=='primary')&(sec['sample']==sample)]
 if ref.empty:continue
 v=v0[(v0.region=='Liver macrometastasis tumor')&(v0.Layer2=='Liver_metastatic carcinoma')].copy();v['near']=v.distance.between(50,200).astype(float);v=v[v.near.eq(1)|v.distance.between(500,1500)].copy();counts=v.groupby('macro_component').near.agg(['sum','count']);good=counts.index[(counts['sum']>=5)&((counts['count']-counts['sum'])>=5)];v=v[v.macro_component.isin(good)].dropna(subset=covs+targets+['neighbor_mural_200'])
 assert len(v)==int(ref.iloc[0]['n']);z=v[covs].to_numpy();active=z.std(0)>1e-10;z=(z[:,active]-z[:,active].mean(0))/z[:,active].std(0);comp=pd.get_dummies(v.macro_component,drop_first=True,dtype=float).to_numpy();a=v.near.to_numpy();nm=v.neighbor_mural_200.to_numpy();nm=(nm-nm.mean())/nm.std();Y=v[targets].to_numpy();base=np.c_[np.ones(len(v)),z,comp]
 for endpoint,X in [('E1',np.c_[base,a]),('E2',np.c_[base,a,nm,a*nm])]:
  beta=lstsq(X,Y,lapack_driver='gelsy')[0][-1]
  r=ref[ref.endpoint==endpoint].set_index('gene').loc[targets];err('primary_QR_coefficients',beta,r.effect);err('primary_SD_coefficients',beta/Y.std(0),r.effect_SD)
  for gene,b in zip(targets,beta):rows.append(dict(sample=sample,endpoint=endpoint,gene=gene,independent_effect=float(b)))
pd.DataFrame(rows).to_csv(O/'results/independent_primary_coefficients.tsv',sep='\t',index=False)
calc=sec.groupby(['variant','endpoint','gene','patient'])[['effect','effect_SD']].mean();saved=pt.set_index(['variant','endpoint','gene','patient']);err('patient_aggregation',calc.loc[saved.index].effect,saved.effect)
def exact_p(x):
 n=len(x);dist=np.array([0.])
 for a in x:dist=np.r_[dist+a,dist-a]
 return np.mean(abs(dist)>=abs(np.sum(x))-1e-11)
for _,r in summary.iterrows():
 x=pt[(pt.variant==r.variant)&(pt.endpoint==r.endpoint)&(pt.gene==r.gene)].effect.to_numpy();assert len(x)==r.n;err('summary_mean',x.mean(),r['mean']);err('summary_P',exact_p(x),r.P)
def bh(p):
 p=np.array(p);order=np.argsort(p);adj=np.empty(len(p));adj[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return adj
for fam,g in summary.groupby('family'):err('BH_summary',bh(g.P.fillna(1)),g.FDR)
for _,r in sp.iterrows():
 q=pt[(pt.variant==r.variant)&(pt.endpoint==r.endpoint)].pivot(index='patient',columns='gene',values='effect_SD').reindex(columns=targets).dropna();delta=q.CGNL1-q[targets[1:]].median(axis=1);err('specificity_mean',delta.mean(),r['mean']);err('specificity_P',exact_p(delta.to_numpy()),r.P)
# Check balanced weighting results directly from recorded weights and independently refit E2.
w=pd.read_csv(O/'results/overlap_weights.tsv.gz',sep='\t');ov=pd.read_csv(O/'results/overlap_section_effects.tsv',sep='\t');balance=[]
for sample,ref in ov.groupby('sample'):
 ww=w[(w['sample']==sample)&(w.status=='qualified')];v=d[d['sample']==sample].merge(ww[['barcode','near','weight']],on='barcode',validate='one_to_one');a=v.near.to_numpy();weight=v.weight.to_numpy();near=a==1;Y=v[targets].to_numpy();b=np.average(Y[near],axis=0,weights=weight[near])-np.average(Y[~near],axis=0,weights=weight[~near]);r=ref[ref.endpoint=='E1'].set_index('gene').loc[targets];err('overlap_E1',b,r.effect)
 Z=v[covs].to_numpy();active=Z.std(0)>1e-10;Z=(Z[:,active]-Z[:,active].mean(0))/Z[:,active].std(0);comp=pd.get_dummies(v.macro_component,drop_first=True,dtype=float).to_numpy();nm=v.neighbor_mural_200.to_numpy();nm=(nm-nm.mean())/nm.std();X=np.c_[np.ones(len(v)),Z,comp,a,nm,a*nm];b=lstsq(X*np.sqrt(weight)[:,None],Y*np.sqrt(weight)[:,None],lapack_driver='gelsy')[0][-1];r=ref[ref.endpoint=='E2'].set_index('gene').loc[targets];err('overlap_E2',b,r.effect)
 delta=np.average(Z[near],axis=0,weights=weight[near])-np.average(Z[~near],axis=0,weights=weight[~near]);balance.append(float(abs(delta).max()))
assert max(balance)<.1
# Confirm original snapshotted project artifacts have not changed.
manifest=json.loads((O/'baseline/SOURCE_MANIFEST.json').read_text());unchanged=all(hashlib.sha256(Path(r['source']).read_bytes()).hexdigest()==r['sha256'] for r in manifest)
assert unchanged
passed=all(v<1e-8 for v in errors.values());assert passed,errors
record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':passed,'errors':errors,'primary_section_gene_coefficients':len(rows),'all_patient_aggregations_checked':len(pt),'all_summary_rows_checked':len(summary),'max_overlap_SMD':max(balance),'original_snapshots_unchanged':unchanged,'frozen_spec_unchanged':True,'scope':'independent primary QR refit, all section-to-patient aggregation and primary/sensitivity summary P/FDR; specificity reconstruction; overlap E1/E2 refit and measured-mean balance; raw-count extraction cross-checked previously against existing audited counts','not_checked':'new independent cohort biological replication, genotype cross-study deduplication, cell-intrinsic origin and manual histopathology reannotation'}
(O/'NUMERICAL_QA.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record,indent=2))
