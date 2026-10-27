from pathlib import Path
import json,gzip,datetime,sys
import numpy as np,pandas as pd
from scipy import sparse,io,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];B=R.parent;D=R/'04_patient_state';C=json.loads((R/'frozen/EXECUTION_CONFIG.json').read_text());old=B/'cohort_verification_20260924';reg=pd.read_csv(old/'qualified_endothelial_patients.tsv',sep='\t');qualified=set(reg.donor)
rows=[];matched=[];inputs=[]
def frozen_endo(X,sy):
 lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;qc=(lib>=500)&(ng>=200)&(mt<=25);scores=[];hits=[];labels=list(C['lineage_markers'])
 for gs in C['lineage_markers'].values():
  z=X[np.isin(sy,gs)].toarray();scores.append(np.log1p(z/np.maximum(lib,1)*1e4).mean(axis=0));hits.append((z>0).sum(axis=0))
 scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2];j=labels.index('endothelial');return qc&(best==j)&(hits[:,j]>=2)&(margin>=.4)
def profile(dataset,donor,sample,X,sy,endo):
 X=X.tocsr();X.eliminate_zeros();assert (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
 if donor not in qualified:return
 k=np.flatnonzero(endo)
 if len(k)<50:return
 lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);assert (lib[k]>0).all();idx=[]
 for g in ['CGNL1']+C['junction_panel']:
  z=np.flatnonzero(sy==g);assert len(z)==1,(dataset,g);idx.append(z[0])
 small=X[idx][:,k].toarray();pb=np.log1p(small.sum(axis=1)/lib[k].sum()*1e4);cell=np.log1p(small/lib[k]*1e4)
 base=dict(dataset=dataset,donor=donor,sample=sample,endo_cells=len(k),CGNL1_positive=int((small[0]>0).sum()),CGNL1_pb=pb[0],junction_pb=float(pb[1:].mean()),mean_logUMI=float(np.log1p(lib[k]).mean()),mean_logGenes=float(np.log1p(ng[k]).mean()),CGNL1_mean_cell=float(cell[0].mean()),junction_mean_cell=float(cell[1:].mean()))
 for j,g in enumerate(['CGNL1']+C['junction_panel']):base[g+'_pb']=pb[j]
 rows.append(base)
 cg=X[idx[0]].toarray().ravel();pos=k[cg[k]>0];neg=k[cg[k]==0];feat=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for i in sorted(pos,key=lambda v:(cg[v],v)):
  can=np.array([v for v in neg if v not in used],int)
  if not len(can):break
  delta=np.abs(feat[can]-feat[i]);ok=(delta<=.25).all(axis=1)
  if ok.any():
   can=can[ok];delta=delta[ok];j=can[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((i,j))
 match=dict(dataset=dataset,donor=donor,sample=sample,pairs=len(pairs),eligible=len(pairs)>=10,effect=np.nan)
 if len(pairs)>=10:
  a,b=np.array(pairs).T;allix=np.r_[a,b];v=np.log1p(X[idx[1:]][:,allix].toarray()/lib[allix]*1e4).mean(axis=0);match['effect']=float((v[:len(a)]-v[len(a):]).mean());match['mean_abs_logUMI_difference']=float(np.abs(feat[a,0]-feat[b,0]).mean())
 matched.append(match);print('PROFILE',donor,sample,len(k),len(pairs),flush=True)
def load_matrix(paths):
 f={k:next(v for v in paths if k in v.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(f['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(f['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
 with gzip.open(f['matrix'],'rb') as h,threadpool_limits(limits=2):X=io.mmread(h).tocsr()
 X.eliminate_zeros();assert X.shape==(len(sy),len(bc));assert len(set(bc))==len(bc)
 inputs.extend(map(str,paths));return X,sy,bc
# discovery, original count matrix and author metadata order from earlier validated conversion
f=B/'deep_v2/results/cnv/counts_gene_by_cell.npz';X=sparse.load_npz(f).tocsr();sy=pd.read_csv(B/'deep_v2/results/cnv/gene_names.tsv',sep='\t').gene.to_numpy(dtype=str);meta=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);assert X.shape==(len(sy),len(meta));endo=frozen_endo(X,sy);inputs.append(str(f))
for donor,g in meta[meta.organs=='LCT'].groupby('patients'):
 ix=meta.index.get_indexer(g.index);profile('GSE225857','GSE225857:'+donor,donor,X[:,ix],sy,endo[ix])
del X
# Previously frozen external cell identities; no reannotation on target expression.
cfg=json.loads((B/'spatial_v3/background_validation/final_exploration_20260924_v2/CONFIG.json').read_text());I=Path(cfg['input_run']);S=Path(cfg['source_run'])
for s in json.loads((I/'SAMPLES.json').read_text()):
 if s['patient_label_provisional'] not in qualified:continue
 if s['dataset']=='GSE231559' and not (s['title'].startswith('L') and not s['title'].endswith('N')):continue
 X,sy,bc=load_matrix([I/f for f in s['files']]);lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;qc=(lib>=500)&(ng>=200)&(mt<=25);X=X[:,qc];bc=bc[qc]
 ann=pd.read_csv(S/'cells'/(s['sample']+'.tsv.gz'),sep='\t').set_index('barcode').reindex(bc);assert ann.label.notna().all()
 profile(s['dataset'],s['patient_label_provisional'],s['sample'],X,sy,ann.label.to_numpy()=='endothelial')
# Additional prior-qualified three-patient cohort, using the same frozen markers.
for s in json.loads((old/'GSE245552_LM_SAMPLES.json').read_text()):
 donor='GSE245552:'+s['donor']
 if donor not in qualified:continue
 X,sy,bc=load_matrix([old/'GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]);profile('GSE245552',donor,s['sample'],X,sy,frozen_endo(X,sy))
q=pd.DataFrame(rows);q.to_csv(D/'lesion_pseudobulk_profiles.tsv',sep='\t',index=False)
p=q.groupby(['dataset','donor'],as_index=False).mean(numeric_only=True);p['qualifying_lesions']=q.groupby(['dataset','donor']).size().to_numpy();assert set(p.donor)==qualified and len(p)==14
p['role']=np.where(p.dataset=='GSE225857','discovery_previously_explored','external_previously_explored');p.to_csv(D/'patient_pseudobulk_profiles.tsv',sep='\t',index=False)
match=pd.DataFrame(matched);match.to_csv(D/'matched_lesion_effects.tsv',sep='\t',index=False);mp=match[match.eligible].groupby(['dataset','donor'],as_index=False).effect.mean();mp.to_csv(D/'matched_patient_effects_descriptive.tsv',sep='\t',index=False)
def fit(v,depth=False):
 x=v.CGNL1_pb.to_numpy();y=v.junction_pb.to_numpy();assert np.std(x,ddof=1)>0 and np.std(y,ddof=1)>0
 design=pd.DataFrame({'CGNL1':(x-x.mean())/x.std(ddof=1)},index=v.index);design=pd.concat([design,pd.get_dummies(v.dataset,prefix='study',drop_first=True,dtype=float)],axis=1)
 if depth:design['logUMI']=(v.mean_logUMI-v.mean_logUMI.mean())/v.mean_logUMI.std(ddof=1)
 design=sm.add_constant(design,has_constant='add');assert np.linalg.matrix_rank(design.to_numpy())==design.shape[1];mod=sm.OLS((y-y.mean())/y.std(ddof=1),design).fit();assert mod.df_resid>=3
 ci=mod.conf_int().loc['CGNL1'];return dict(n_patients=len(v),n_studies=v.dataset.nunique(),beta=float(mod.params['CGNL1']),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(mod.pvalues['CGNL1']),residual_df=float(mod.df_resid),design_condition=float(np.linalg.cond(design.to_numpy())),R2=float(mod.rsquared))
tests=[]
for name,v in [('all_previously_explored',p),('external_previously_explored',p[p.dataset!='GSE225857'])]:tests.append(dict(model=name,**fit(v),family='A05_continuous_primary_2',scope='exploratory_association_not_new_confirmation'))
test=pd.DataFrame(tests);test['FDR']=multipletests(test.P,method='fdr_bh')[1];test['nominal']=test.P<.05;test['FDR_significant']=test.FDR<.05;test.to_csv(D/'A05_primary_models.tsv',sep='\t',index=False)
sens=pd.DataFrame([dict(model='all_plus_depth',**fit(p,True),family='prespecified_depth_sensitivity_single')]);sens.to_csv(D/'A05_depth_sensitivity.tsv',sep='\t',index=False)
loo=[]
for donor in p.donor:
 try:loo.append(dict(omitted_patient=donor,**fit(p[p.donor!=donor])))
 except Exception as e:loo.append(dict(omitted_patient=donor,status='not_estimable',reason=str(e)))
pd.DataFrame(loo).to_csv(D/'leave_one_patient_out.tsv',sep='\t',index=False)
by=[]
for ds,v in p.groupby('dataset'):
 rho=stats.spearmanr(v.CGNL1_pb,v.junction_pb).statistic if len(v)>1 else np.nan
 pv=stats.spearmanr(v.CGNL1_pb,v.junction_pb).pvalue if len(v)>=4 else np.nan
 by.append(dict(dataset=ds,n=len(v),spearman_rho=rho,P=pv,family='A05_per_study_estimable',evidence='exploratory' if len(v)>=4 else 'descriptive'))
by=pd.DataFrame(by);by['FDR']=np.nan;mask=by.P.notna()
if mask.any():by.loc[mask,'FDR']=multipletests(by.loc[mask,'P'],method='fdr_bh')[1]
by.to_csv(D/'per_study_association.tsv',sep='\t',index=False)
(D/'INPUT_PATHS.json').write_text(json.dumps(sorted(set(inputs)),indent=2));summary={'exit_code':0,'scope':'A05 initial model; A04 subtype and A06 control-gene sensitivity remain pending','patient_n':len(p),'discovery_n':int((p.dataset=='GSE225857').sum()),'external_previously_explored_n':int((p.dataset!='GSE225857').sum()),'matched_eligible_patients':len(mp),'primary_results':test.to_dict('records'),'finished':datetime.datetime.now().astimezone().isoformat()};(D/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
