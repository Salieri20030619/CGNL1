import pathlib,json,numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_routes';plan=json.loads((O/'PLAN.json').read_text());gs=plan['genes']
m=pd.read_csv(V/'03_processed/GSE234129_annotation_qc.tsv.gz',sep='\t');bar=pd.read_csv(B/'02_raw/GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0];ft=pd.read_csv(B/'02_raw/GSE234129_features.tsv.gz',sep='\t',header=None);sy=ft.iloc[:,1 if ft.shape[1]>1 else 0].astype(str).values
with threadpool_limits(limits=2):X=io.mmread(B/'02_raw/GSE234129_count_matrix.mtx.gz').tocsr()
ix=pd.Index(bar).get_indexer(m.cell_barcodes);assert (ix>=0).all();X=X[:,ix];umi=np.asarray(X.sum(axis=0)).ravel();gi=[np.flatnonzero(sy==g) for g in gs];assert all(len(z)==1 for z in gi);Y=X[np.array(gi).ravel()].toarray().astype(float);N=Y/np.maximum(umi,1)[None,:]*1e4
m['umi_recomputed']=umi;m.to_csv(O/'source_cell_metadata.tsv.gz',sep='\t',index=False)
# Direct sample and parent-group expression/capture summaries.
rows=[]
for level in ['broad_type','celltype']:
 for (pt,sample,ct),ix in m.groupby(['patient','sample',level]).groups.items():
  ix=np.asarray(list(ix));site=sample.split('-')[-1];n_sample=sum(m['sample']==sample)
  for j,g in enumerate(gs):rows.append({'level':level,'patient':pt,'sample':sample,'site':site,'celltype':ct,'parent':m.loc[ix,'broad_type'].iloc[0],'gene':g,'n_cells':len(ix),'sample_cells':n_sample,'fraction_of_captured_cells':len(ix)/n_sample,'counts':Y[j,ix].sum(),'library_size':umi[ix].sum(),'mean_normalized_1e4':N[j,ix].mean(),'CPM':Y[j,ix].sum()/max(umi[ix].sum(),1)*1e6,'detection_fraction':(Y[j,ix]>0).mean(),'median_umi':np.median(umi[ix])})
a=pd.DataFrame(rows);a['log2CPM']=np.log2(a.CPM+1);a.to_csv(O/'paired_gene_expression_by_subtype.tsv',sep='\t',index=False)
pt_sites=pd.crosstab(m.patient,m.site);pts=pt_sites.index[(pt_sites.Ca>0)&(pt_sites.As>0)].tolist();assert len(pts)==4
pairs=[];props=[];decomp=[];detail=[]
for pt in pts:
 for level in ['broad_type','celltype']:
  z=m[(m.patient==pt)&m.site.isin(['Ca','As'])];ct=pd.crosstab(z[level],z.site).reindex(index=sorted(m[m.patient.isin(pts)&m.site.isin(['Ca','As'])][level].unique()),columns=['Ca','As'],fill_value=0)
  for c in ct.index:props.append({'patient':pt,'level':level,'celltype':c,'primary_n':ct.loc[c,'Ca'],'ascites_n':ct.loc[c,'As'],'primary_fraction':ct.loc[c,'Ca']/ct.Ca.sum(),'ascites_fraction':ct.loc[c,'As']/ct.As.sum(),'delta_fraction':ct.loc[c,'As']/ct.As.sum()-ct.loc[c,'Ca']/ct.Ca.sum()})
 for level in ['broad_type','celltype']:
  for cut in [20,10]:
   z=a[(a.patient==pt)&a.site.isin(['Ca','As'])&(a.level==level)&(a.n_cells>=cut)]
   for (g,c),v in z.groupby(['gene','celltype']):
    if set(v.site)!=set(['Ca','As']):continue
    p=v.set_index('site').loc['Ca'];q=v.set_index('site').loc['As'];pairs.append({'patient':pt,'level':level,'min_cells':cut,'gene':g,'celltype':c,'primary_n':p.n_cells,'ascites_n':q.n_cells,'primary_log2CPM':p.log2CPM,'ascites_log2CPM':q.log2CPM,'delta_log2CPM':q.log2CPM-p.log2CPM,'primary_mean_norm':p.mean_normalized_1e4,'ascites_mean_norm':q.mean_normalized_1e4,'delta_mean_norm':q.mean_normalized_1e4-p.mean_normalized_1e4,'delta_detection_fraction':q.detection_fraction-p.detection_fraction})
 # Whole captured population and within two common lineages.
 for parent in ['All_captured','Macrophage','Lymphoid']:
  mm=m[(m.patient==pt)&m.site.isin(['Ca','As'])];mm=mm if parent=='All_captured' else mm[mm.broad_type==parent]
  level='broad_type' if parent=='All_captured' else 'celltype'
  for cut in [20,10]:
   groups={site:{ct:np.asarray(list(ix)) for ct,ix in mm[mm.site==site].groupby(level).groups.items()} for site in ['Ca','As']};tot={s:sum(len(ix) for ix in groups[s].values()) for s in groups};eligible=[c for c in set(groups['Ca'])&set(groups['As']) if min(len(groups[s][c]) for s in groups)>=cut]
   for j,g in enumerate(gs):
    mu={s:N[j,mm[mm.site==s].index].mean() for s in groups};state=0.;composition=0.
    for c in sorted(eligible):
     ip,ia=groups['Ca'][c],groups['As'][c];pp,pa=len(ip)/tot['Ca'],len(ia)/tot['As'];ep,ea=N[j,ip].mean(),N[j,ia].mean();ds=(pa+pp)/2*(ea-ep);dc=(ea+ep)/2*(pa-pp);assert abs((pa*ea-pp*ep)-(ds+dc))<1e-10;state+=ds;composition+=dc
     detail.append({'patient':pt,'parent':parent,'min_cells':cut,'gene':g,'celltype':c,'primary_fraction':pp,'ascites_fraction':pa,'primary_mean_norm':ep,'ascites_mean_norm':ea,'state_component':ds,'composition_component':dc})
    delta=mu['As']-mu['Ca'];residual=delta-state-composition;decomp.append({'patient':pt,'parent':parent,'min_cells':cut,'gene':g,'primary_mean_norm':mu['Ca'],'ascites_mean_norm':mu['As'],'total_delta':delta,'shared_state':state,'shared_composition':composition,'unresolved_rare_or_absent':residual,'shared_types':';'.join(sorted(eligible)),'primary_coverage':sum(len(groups['Ca'][c]) for c in eligible)/tot['Ca'],'ascites_coverage':sum(len(groups['As'][c]) for c in eligible)/tot['As']})
p=pd.DataFrame(pairs);p.to_csv(O/'within_type_paired_effects.tsv',sep='\t',index=False);summary=p.groupby(['level','min_cells','gene','celltype']).agg(n_pairs=('patient','nunique'),median_log2CPM_delta=('delta_log2CPM','median'),positive=('delta_log2CPM',lambda x:int((x>0).sum())),negative=('delta_log2CPM',lambda x:int((x<0).sum())),zero=('delta_log2CPM',lambda x:int((x==0).sum())),median_mean_norm_delta=('delta_mean_norm','median'),mean_norm_positive=('delta_mean_norm',lambda x:int((x>0).sum())),median_detection_delta=('delta_detection_fraction','median')).reset_index();summary.to_csv(O/'within_type_summary.tsv',sep='\t',index=False)
pr=pd.DataFrame(props);pr.to_csv(O/'captured_cell_proportion_paired.tsv',sep='\t',index=False);pq=pr.groupby(['level','celltype']).agg(n_pairs=('patient','nunique'),primary_fraction_median=('primary_fraction','median'),ascites_fraction_median=('ascites_fraction','median'),median_delta_fraction=('delta_fraction','median'),positive=('delta_fraction',lambda x:int((x>0).sum())),negative=('delta_fraction',lambda x:int((x<0).sum()))).reset_index();pq.to_csv(O/'captured_cell_proportion_summary.tsv',sep='\t',index=False)
d=pd.DataFrame(decomp);d.to_csv(O/'state_composition_decomposition.tsv',sep='\t',index=False);pd.DataFrame(detail).to_csv(O/'decomposition_by_shared_type.tsv',sep='\t',index=False);assert np.max(abs(d.total_delta-d.shared_state-d.shared_composition-d.unresolved_rare_or_absent))<1e-10
(O/'PAIRED_COMPLETE.json').write_text(json.dumps({'exit_code':0,'patients':pts,'cells':len(m),'decomposition_identity_max_error':float(np.max(abs(d.total_delta-d.shared_state-d.shared_composition-d.unresolved_rare_or_absent))),'n_inferential_tests':0},indent=2));print(summary[(summary.min_cells==20)&((summary.level=='broad_type')|summary.celltype.str.startswith('TAM'))].to_string(index=False));print(pq[(pq.level=='broad_type')|pq.celltype.str.startswith('TAM')].to_string(index=False));print(d[(d.min_cells==20)&(d.parent=='Macrophage')].to_string(index=False))

