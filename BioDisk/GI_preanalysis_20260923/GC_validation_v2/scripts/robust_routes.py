import pathlib,json,numpy as np,pandas as pd
from scipy import stats
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2/candidate_routes');plan=json.loads((O/'PLAN.json').read_text());(O/'ROBUSTNESS_PLAN.json').write_text(json.dumps({'trigger':'Two phagolysosomal associations passed primary and covariate sensitivity FDR; investigate single-patient/marker dependence and low-cell-count compartments before interpreting.','analyses':['leave-one-patient descriptive correlation ranges','six leave-one-program-marker-out models per two observed associations, BH across12 sensitivity tests','min_cells>=50 repeats full54 family, raw and adjusted each within54','macrophage-subtype captured fractions with all four patients incl zeros']},indent=2))
a=pd.read_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t');p=pd.read_csv(O/'external_program_scores.tsv',sep='\t');r0=pd.read_csv(O/'external_state_program_correlations.tsv',sep='\t')
def partial(z):
 n=len(z);dum=pd.get_dummies(z.sample_type,drop_first=True,dtype=float).values;C=np.column_stack([np.ones(n),dum,stats.rankdata(np.log1p(z.n_cells)),stats.rankdata(np.log1p(z.median_umi))]);rk=np.linalg.matrix_rank(C);x=stats.rankdata(z.log2CPM);y=stats.rankdata(z.score);x-=C@np.linalg.lstsq(C,x,rcond=None)[0];y-=C@np.linalg.lstsq(C,y,rcond=None)[0];df=n-rk-1
 if min(np.std(x),np.std(y))<1e-10:return np.nan,np.nan
 rr=np.corrcoef(x,y)[0,1];pv=2*stats.t.sf(abs(rr)*np.sqrt(df/max(1-rr*rr,1e-15)),df);return rr,pv
sel=r0[(r0.FDR_54<.05)&(r0.FDR_partial_54<.05)]
loo=[];log=[]
for _,r in sel.iterrows():
 ct,g,prog=r.celltype,r.gene,r.program;z=a[(a.celltype==ct)&(a.gene==g)&(a.n_cells>=20)][['patient','sample_type','n_cells','median_umi','log2CPM']].merge(p[(p.celltype==ct)&(p.program==prog)][['patient','score']],on='patient',validate='one_to_one')
 for pt in z.patient:
  zz=z[z.patient!=pt];rr=stats.spearmanr(zz.log2CPM,zz.score).statistic;pr,_=partial(zz);loo.append({'gene':g,'celltype':ct,'program':prog,'omitted_patient':pt,'rho':rr,'partial_rho':pr})
 for omit in plan['programs'][prog]:
  rest=[q for q in plan['programs'][prog] if q!=omit];s=a[(a.celltype==ct)&a.gene.isin(rest)].groupby('patient').log2CPM.mean().rename('score');zz=z.drop(columns='score').merge(s,on='patient');co=stats.spearmanr(zz.log2CPM,zz.score);pr,pv=partial(zz);log.append({'gene':g,'celltype':ct,'program':prog,'omitted_marker':omit,'n':len(zz),'rho':co.statistic,'P':co.pvalue,'partial_rho':pr,'P_partial':pv})
l=pd.DataFrame(log)
for pcol,qcol in [('P','FDR_12_sensitivity'),('P_partial','FDR_partial_12_sensitivity')]:l[qcol]=stats.false_discovery_control(l[pcol].fillna(1).values,method='bh')
l.to_csv(O/'leave_one_program_gene_out.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'leave_one_patient_correlations.tsv',sep='\t',index=False)
# Full family at larger minimum cell count.
rows=[]
for _,r in r0.iterrows():
 ct,g,prog=r.celltype,r.gene,r.program;z=a[(a.celltype==ct)&(a.gene==g)&(a.n_cells>=50)][['patient','sample_type','n_cells','median_umi','log2CPM']].merge(p[(p.celltype==ct)&(p.program==prog)][['patient','score']],on='patient',validate='one_to_one');n=len(z);v={'gene':g,'celltype':ct,'program':prog,'n':n,'rho':np.nan,'P':np.nan,'partial_rho':np.nan,'P_partial':np.nan}
 if n>=10 and z.log2CPM.nunique()>2 and z.score.nunique()>2:
  c=stats.spearmanr(z.log2CPM,z.score);pr,pv=partial(z);v.update(rho=c.statistic,P=c.pvalue,partial_rho=pr,P_partial=pv)
 rows.append(v)
r=pd.DataFrame(rows)
for pcol,qcol in [('P','FDR_54_min50'),('P_partial','FDR_partial_54_min50')]:r[qcol]=stats.false_discovery_control(r[pcol].fillna(1).values,method='bh');r.loc[r[pcol].isna(),qcol]=np.nan
r.to_csv(O/'program_correlations_min50.tsv',sep='\t',index=False)
m=pd.read_csv(O/'source_cell_metadata.tsv.gz',sep='\t');m=m[(m.broad_type=='Macrophage')&m.site.isin(['Ca','As'])];pts=sorted(set(m[m.site=='Ca'].patient)&set(m[m.site=='As'].patient));sub=sorted(m.celltype.unique());rows=[]
for pt in pts:
 z=m[m.patient==pt];ct=pd.crosstab(z.celltype,z.site).reindex(index=sub,columns=['Ca','As'],fill_value=0)
 for c in sub:rows.append({'patient':pt,'subtype':c,'primary_n':ct.loc[c,'Ca'],'ascites_n':ct.loc[c,'As'],'primary_fraction_in_macrophages':ct.loc[c,'Ca']/ct.Ca.sum(),'ascites_fraction_in_macrophages':ct.loc[c,'As']/ct.As.sum(),'delta_fraction':ct.loc[c,'As']/ct.As.sum()-ct.loc[c,'Ca']/ct.Ca.sum()})
f=pd.DataFrame(rows);f.to_csv(O/'macrophage_subtype_fractions.tsv',sep='\t',index=False);q=f.groupby('subtype').agg(n=('patient','nunique'),primary_median=('primary_fraction_in_macrophages','median'),ascites_median=('ascites_fraction_in_macrophages','median'),median_delta=('delta_fraction','median'),positive=('delta_fraction',lambda x:int((x>0).sum())),negative=('delta_fraction',lambda x:int((x<0).sum()))).reset_index();q.to_csv(O/'macrophage_subtype_fraction_summary.tsv',sep='\t',index=False)
# External abundance and gene transcript shares: absent categories are true zero in captured material only.
src=pd.read_csv(O/'external_candidate_source_shares.tsv',sep='\t');types=sorted(src.celltype.unique());patients=sorted(src.patient.unique());filled=[]
for g in plan['genes']:
 for pt in patients:
  z=src[(src.gene==g)&(src.patient==pt)].set_index('celltype')
  for ct in types:
   if ct in z.index:rr=z.loc[ct];filled.append({'gene':g,'patient':pt,'celltype':ct,'cell_fraction':rr.cell_fraction,'gene_count_share':rr.gene_count_share,'n_cells':rr.n_cells})
   else:filled.append({'gene':g,'patient':pt,'celltype':ct,'cell_fraction':0.,'gene_count_share':0.,'n_cells':0})
ff=pd.DataFrame(filled);ff.to_csv(O/'external_source_shares_zero_filled.tsv',sep='\t',index=False);qq=ff.groupby(['gene','celltype']).agg(n_patients=('patient','nunique'),median_cell_fraction=('cell_fraction','median'),median_gene_count_share=('gene_count_share','median'),patients_ge20=('n_cells',lambda x:int((x>=20).sum()))).reset_index();qq.to_csv(O/'external_source_shares_zero_filled_summary.tsv',sep='\t',index=False)
print(q.to_string(index=False));print(l.to_string(index=False));print(r[(r.gene.isin(['CDH23','COPZ2']))&(r.celltype=='Macrophage_enriched')&(r.program=='phagolysosomal')].to_string(index=False));print(pd.DataFrame(loo).groupby('gene')[['rho','partial_rho']].agg(['min','max']).to_string());(O/'ROBUSTNESS_COMPLETE.json').write_text(json.dumps({'exit_code':0,'selected_associations':len(sel),'leave_one_marker_models':len(l),'min50_family':len(r),'macrophage_paired_patients':len(pts)},indent=2))
