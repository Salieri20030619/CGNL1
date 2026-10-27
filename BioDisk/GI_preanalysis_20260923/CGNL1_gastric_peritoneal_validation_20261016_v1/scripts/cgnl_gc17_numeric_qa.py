from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import json,datetime,hashlib,itertools
import numpy as np,pandas as pd
from scipy import stats
O=Path.cwd();errors={}
def err(k,a,b):errors[k]=max(errors.get(k,0),float(np.max(np.abs(np.array(a)-np.array(b)))))
def hc3(y,X,idx):
 beta=np.linalg.lstsq(X,y,rcond=None)[0];inv=np.linalg.inv(X.T@X);lev=np.sum((X@inv)*X,axis=1);res=y-X@beta;cov=inv@(X.T@((res/(1-lev))[:,None]**2*X))@inv;se=np.sqrt(cov[idx,idx]);df=len(y)-X.shape[1];P=2*stats.t.sf(abs(beta[idx]/se),df);return beta[idx],P
d=pd.read_csv(O/'results/bulk_patient_deltas.tsv',sep='\t');models=pd.read_csv(O/'results/bulk_models.tsv',sep='\t');panels=list(json.loads((O/'IMPLEMENTATION_FREEZE.json').read_text())['bulk_panels'])
for name,covs in [('M0',[]),('M1',['preservation_delta']),('M2',['preservation_delta']+panels)]:
 X=np.column_stack([np.ones(len(d))]+[d[c] for c in covs]);b,P=hc3(d.CGNL1.to_numpy(),X,0);r=models[models.model==name].iloc[0];err('bulk_HC3_effect',b,r.effect);err('bulk_HC3_P',P,r.P)
bres=json.loads((O/'results/BULK_RESULT.json').read_text());err('paired_signed_rank',stats.wilcoxon(d.CGNL1,method='exact').pvalue,bres['P_exact_signed_rank'])
m=pd.read_csv(O/'processed/original31_sample_values.tsv',sep='\t');m=m[m.group_from_filename!='healthy_unpaired'];models=pd.read_csv(O/'results/original31_models.tsv',sep='\t');pc=m.group_from_filename.eq('with_PC').astype(float).to_numpy();ver=pd.get_dummies(m.program_version,drop_first=True,dtype=float)
for name,covs in [('M0',[]),('M1',list(ver.columns)),('M2',list(ver.columns)+panels)]:
 X=np.column_stack([np.ones(len(m)),pc]+[ver[c] if c in ver else m[c] for c in covs]);b,P=hc3(m.CGNL1.to_numpy(),X,1);r=models[models.model==name].iloc[0];err('original31_HC3_effect',b,r.effect);err('original31_HC3_P',P,r.P)
orig=json.loads((O/'results/ORIGINAL31_RESULT.json').read_text());a=m[m.group_from_filename=='with_PC'].CGNL1;b=m[m.group_from_filename=='without_PC'].CGNL1;se=np.sqrt(a.var()/len(a)+b.var()/len(b));df=(a.var()/len(a)+b.var()/len(b))**2/((a.var()/len(a))**2/(len(a)-1)+(b.var()/len(b))**2/(len(b)-1));P=2*stats.t.sf(abs((a.mean()-b.mean())/se),df);err('original31_Welch_P',P,orig['P_Welch'])
def bh(x):
 x=np.array(x);order=np.argsort(x);out=np.empty(len(x));out[order]=np.minimum(1,np.minimum.accumulate((x[order]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]);return out
for fname in ['bulk_models.tsv','original31_models.tsv']:
 tab=pd.read_csv(O/'results'/fname,sep='\t');err('BH_'+fname,bh(tab.P),tab.FDR_three_models)
prof=pd.read_csv(O/'results/identity_source_profiles.tsv',sep='\t');contr=pd.read_csv(O/'results/identity_source_contrasts.tsv',sep='\t')
for _,r in contr.iterrows():
 z=prof[prof.n_cells>=r.min_cells];other=r.contrast.split('-minus-')[1];p=z[z.group.isin(['strict_mesothelial',other])].pivot(index='patient',columns='group',values='CGNL1').dropna();delta=p.strict_mesothelial-p[other];assert len(p)==r.n;err('source_effect',delta.mean(),r['mean'])
 if np.isfinite(r.P):
  v=delta[delta!=0].to_numpy();ranks=stats.rankdata(abs(v));ss=np.array([0.])
  for k in ranks:ss=np.r_[ss,ss+k]
  obs=min(ranks[v>0].sum(),ranks[v<0].sum());P=min(1,2*np.mean(ss<=obs+1e-10));err('source_P',P,r.P)
for cutoff,tab in contr.groupby('min_cells'):
 tab=tab.dropna(subset=['P'])
 if len(tab):err('source_BH',bh(tab.P),tab.FDR)
good=prof.n_cells>0;expected=np.log1p(prof.loc[good,'CGNL1_count']/prof.loc[good,'total_UMI']*1e4);err('profile_normalization',expected,prof.loc[good,'CGNL1']);assert ((prof.CGNL1_detected<=prof.n_cells)&(prof.CGNL1_count>=prof.CGNL1_detected)).all()
states=pd.read_csv(O/'results/mesothelial_state_models.tsv',sep='\t')
for _,r in states[(states.model=='HC3_depth_identity')&(states.status=='estimated')].iterrows():
 z=prof[(prof.group=='strict_mesothelial')&(prof.n_cells>=r.min_cells)];X=np.c_[np.ones(len(z)),z.CGNL1,z.mean_logUMI,z.mesothelial_identity];b,P=hc3(z[r.panel].to_numpy(),X,1);err('state_HC3_effect',b,r.effect);err('state_HC3_P',P,r.P)
assert max(errors.values())<1e-8,errors
receipt=json.loads((O/'metadata/ORIGINAL31_DOWNLOAD_RECEIPTS.json').read_text());assert all(hashlib.sha256((O/'raw/original31'/r['filename']).read_bytes()).hexdigest()==r['expected_sha256'] for r in receipt)
freeze=json.loads((O/'SPEC_FREEZE.json').read_text());assert hashlib.sha256((O/'ANALYSIS_SPECIFICATION.json').read_bytes()).hexdigest()==freeze['sha256']
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'max_errors':errors,'author_CHPSHA_verified':54,'original31_binary_decoders_agree':json.loads((O/'metadata/CHP_DECODE_QUALIFICATION.json').read_text())['parallel_binary_decoders_agree'],'new_expected_counts_match_existing':json.loads((O/'metadata/NEW_BULK_DOWNLOAD_VERIFICATION.json').read_text())['table_identical'],'freeze_unchanged':True,'scope':'independentHC3matrixalgebra/Welch/signedrank/BH;sourceprofilenormalization;rawCHPchecksums;notclinicalindependentlabelverification'};(O/'NUMERICAL_QA.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
# Collapse original author subtypes only for descriptive cell-source summary; preserve original table.
source=pd.read_csv(O/'results/GSE234129_author_source.tsv',sep='\t')
def broad(x):
 for k in ['Fibroblast','Endothelial','Mesothelial','VSMC']:
  if x.startswith(k):return k
 return 'Immune'
source['broad_type']=source.celltype.map(broad);g=source.groupby(['patient','sample','broad_type'],as_index=False)[['n_cells','target_counts','total_UMI','detected']].sum();g['ln1pCP10k']=np.log1p(g.target_counts/g.total_UMI*1e4);g['fraction_detected']=g.detected/g.n_cells;g.to_csv(O/'results/GSE234129_broad_source.tsv',sep='\t',index=False)
print('234_BROAD',g.groupby('broad_type')[['n_cells','target_counts']].sum().to_string());print('MESO_SUMMARY',pd.read_csv(O/'results/identity_source_summary.tsv',sep='\t').to_string(index=False))
