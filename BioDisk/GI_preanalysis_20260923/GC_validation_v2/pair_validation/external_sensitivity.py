import pathlib,json,numpy as np,pandas as pd
from scipy import stats
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2/pair_validation');plan=json.loads((O/'PLAN.json').read_text());a=pd.read_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t');scores=pd.read_csv(O/'external_program_scores.tsv',sep='\t');r0=pd.read_csv(O/'external_state_program_correlations.tsv',sep='\t')
# Sensitivities fixed by PLAN: add gastric differentiation and raise cell-count threshold.
def partial(z,gastric=False):
 n=len(z);dum=pd.get_dummies(z.sample_type,drop_first=True,dtype=float).values;C=np.column_stack([np.ones(n),dum,stats.rankdata(np.log1p(z.n_cells)),stats.rankdata(np.log1p(z.median_umi))]+([stats.rankdata(z.gastric_score)] if gastric else []));rk=np.linalg.matrix_rank(C);x=stats.rankdata(z.log2CPM);y=stats.rankdata(z.score);x=x-C@np.linalg.lstsq(C,x,rcond=None)[0];y=y-C@np.linalg.lstsq(C,y,rcond=None)[0];df=n-rk-1
 if n<10 or df<3 or min(np.std(x),np.std(y))<1e-10:return np.nan,np.nan,df
 rr=float(np.corrcoef(x,y)[0,1]);return rr,float(2*stats.t.sf(abs(rr)*np.sqrt(df/max(1-rr*rr,1e-15)),df)),df
res=[];loo=[];lomo=[]
for cutoff in [20,50]:
 for _,r in r0.iterrows():
  ct,g,prog=r.celltype,r.gene,r.program;x=a[(a.celltype==ct)&(a.gene==g)&(a.n_cells>=cutoff)][['patient','sample_type','n_cells','median_umi','log2CPM']];z=x.merge(scores[(scores.celltype==ct)&(scores.program==prog)][['patient','score']],on='patient',validate='one_to_one');gs=a[(a.celltype==ct)&a.gene.isin(plan['gastric_proxy'])].groupby('patient').log2CPM.mean().rename('gastric_score');z=z.merge(gs,on='patient',validate='one_to_one');n=len(z)
  raw=stats.spearmanr(z.log2CPM,z.score) if n>=10 and z.log2CPM.nunique()>2 and z.score.nunique()>2 else None;pr,pp,df=partial(z);gr,gp,gdf=partial(z,True);res.append({'gene':g,'celltype':ct,'program':prog,'min_cells':cutoff,'n':n,'rho':raw.statistic if raw else np.nan,'P':raw.pvalue if raw else np.nan,'partial_rho':pr,'P_partial':pp,'df':df,'gastric_partial_rho':gr,'P_gastric':gp,'gastric_df':gdf})
  if cutoff==20 and g=='OPLAH' and ct=='Epithelial_like' and prog=='glutathione_cycle':
   for pt in z.patient:
    zz=z[z.patient!=pt];rr=stats.spearmanr(zz.log2CPM,zz.score);pr,_,_=partial(zz);loo.append({'omitted_patient':pt,'rho':rr.statistic,'partial_rho':pr})
   for omit in plan['programs'][prog]:
    rest=[q for q in plan['programs'][prog] if q!=omit];ss=a[(a.celltype==ct)&a.gene.isin(rest)].groupby('patient').log2CPM.mean().rename('score');zz=z.drop(columns='score').merge(ss,on='patient');rr=stats.spearmanr(zz.log2CPM,zz.score);pr,pp,_=partial(zz);lomo.append({'omitted_marker':omit,'rho':rr.statistic,'P':rr.pvalue,'partial_rho':pr,'P_partial':pp})
r=pd.DataFrame(res)
for cutoff in [20,50]:
 ix=r.min_cells==cutoff;assert ix.sum()==24
 for pcol,qcol in [('P','FDR24'),('P_partial','FDR_partial24'),('P_gastric','FDR_gastric24')]:
  r.loc[ix,qcol]=stats.false_discovery_control(r.loc[ix,pcol].fillna(1).values);r.loc[ix&r[pcol].isna(),qcol]=np.nan
r.to_csv(O/'external_sensitivities.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'OPLAH_glutathione_LOO_patients.tsv',sep='\t',index=False);pd.DataFrame(lomo).to_csv(O/'OPLAH_glutathione_LOO_markers_descriptive.tsv',sep='\t',index=False)
# Marker and patient leave-one-out analyses are descriptive sensitivity ranges, not replacement primary tests.
(O/'SENSITIVITY_COMPLETE.json').write_text(json.dumps({'exit_code':0,'family_per_cutoff':24,'rows':len(r),'LOO_patients':len(loo),'LOO_markers':len(lomo),'selection_for_LOO':'OPLAH/epithelial/glutathione was sole raw FDR positive association; technical-adjusted FDR negative retained'},indent=2));print(r[(r.celltype=='Epithelial_like')&(((r.gene=='OPLAH')&(r.program=='glutathione_cycle'))|((r.gene=='RAP1GAP2')&(r.program=='adhesion')))].to_string(index=False));print(pd.DataFrame(loo).agg(['min','max']).to_string());print(pd.DataFrame(lomo).to_string(index=False))
