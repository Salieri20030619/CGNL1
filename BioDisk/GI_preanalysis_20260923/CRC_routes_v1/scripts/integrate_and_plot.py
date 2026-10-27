import pathlib,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');A=O/'results/route1';B=O/'results/route2';F=O/'results/figures';F.mkdir(exist_ok=True)
# Cross-dataset cell source, using all five matched colon patients as a reference, not independent metastatic replication.
g=pd.read_csv(A/'group_metadata.tsv',sep='\t');n=np.load(A/'aggregates.npz');sy=n['genes'];L=n['log2CPM'];screen=pd.read_csv(B/'descriptive_screen.tsv',sep='\t');ss=screen[screen.contrast=='tumor_front_vs_other_tumor'];rows=[]
for gene in sorted(set(ss.gene)):
 gi=np.flatnonzero(sy==gene)
 if len(gi)!=1:continue
 for pt in g.patient.unique():
  vals={}
  for ct in ['Tumor_author','Fibroblast','Endothelial']:
   ix=g[(g.patient==pt)&(g.celltype==ct)&(g.n_cells>=20)].index.values;vals[ct]=float(L[gi[0],ix].mean()) if len(ix) else np.nan
  rows.append({'gene':gene,'patient':pt,**vals,'fibroblast_minus_tumor':vals['Fibroblast']-vals['Tumor_author']})
r=pd.DataFrame(rows);r.to_csv(B/'candidate_cell_source_by_patient.tsv',sep='\t',index=False);src=r.groupby('gene').agg(reference_patients=('patient','nunique'),median_fibroblast_minus_tumor=('fibroblast_minus_tumor','median'),fibroblast_higher=('fibroblast_minus_tumor',lambda x:int((x>0).sum())),median_tumor_log2CPM=('Tumor_author','median'),median_fibroblast_log2CPM=('Fibroblast','median')).reset_index();out=ss.merge(src,on='gene',how='left');out.to_csv(B/'front_screen_with_cell_source.tsv',sep='\t',index=False)
# Does the invasion proxy change co-vary with the non-overlapping stromal proxy across patients? Descriptive rho only n=5.
p=pd.read_csv(B/'program_patient_effects.tsv',sep='\t');w=p[p.contrast=='tumor_front_vs_other_tumor'].pivot(index='patient',columns='program',values='delta');corr=w[['invasion_ECM','stromal_proxy','epithelial_available3']].corr(method='spearman');corr.to_csv(B/'descriptive_program_change_correlations.tsv',sep='\t')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
p1=pd.read_csv(A/'program_patient_effects.tsv',sep='\t');p1=p1[p1.celltype=='Tumor_strict_QC'];he=pd.read_csv(A/'hepatic_proxy_patient_effects.tsv',sep='\t')
for ax,data,names,title in [(axs[0],p1,['invasion_ECM','intestinal_differentiation','proliferation'],'Route 1: tumor-label QC sensitivity'),(axs[1],he,['Tumor_strict_QC','Fibroblast','Endothelial'],'Liver-associated RNA by cell compartment')]:
 col='program' if 'program' in data else 'celltype'
 for i,name in enumerate(names):
  z=data[data[col]==name].sort_values('patient');ax.scatter(np.arange(len(z))*0.035+i-.07,z.delta,s=40,c='#356DA0');ax.hlines(z.delta.median(),i-.18,i+.18,color='#B85532',lw=2)
 ax.axhline(0,color='grey',lw=.8);ax.set_xticks(range(len(names)),[x.replace('_','\n') for x in names]);ax.set_ylabel('Liver minus primary: mean log2(CPM + 1) score');ax.set_title(title,fontsize=11);ax.grid(axis='y',alpha=.15)
fig.suptitle('Five confirmed colon-cancer patients; paired descriptive analysis\nAll received preoperative chemotherapy; no inferential P/FDR',fontsize=12);fig.savefig(F/'route1_patient_effects.png',dpi=180);fig.savefig(F/'route1_patient_effects.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
for ax,ct,title in [(axs[0],'tumor_front_vs_other_tumor','Front vs other tumor regions (5 patients)'),(axs[1],'tumor_front_vs_deep_tumor','Front vs deep tumor regions (4 patients)')]:
 z=p[p.contrast==ct]
 for i,name in enumerate(['invasion_ECM','intestinal_differentiation','stromal_proxy','epithelial_available3']):
  a=z[z.program==name].sort_values('patient');ax.scatter(np.arange(len(a))*.035+i-.07,a.delta,s=40,c='#356DA0');ax.hlines(a.delta.median(),i-.18,i+.18,color='#B85532',lw=2)
 ax.axhline(0,color='grey',lw=.8);ax.set_xticks(range(4),['Invasion / ECM','Intestinal\ndifferentiation','Stromal proxy','Epithelial\n3-gene proxy'],rotation=15);ax.set_ylabel('Patient mean of section-level score differences');ax.set_title(title,fontsize=11);ax.grid(axis='y',alpha=.15)
fig.suptitle('Route 2: pathology-annotated invasive front\nExploratory CRC spatial cohort; individual colon/rectum mapping unresolved',fontsize=12);fig.savefig(F/'route2_patient_effects.png',dpi=180);fig.savefig(F/'route2_patient_effects.svg');plt.close(fig)
a=pd.read_csv(B/'spot_annotation_coordinates.tsv.gz',sep='\t');a=a[a['sample']=='M-ST-29'];mapping={'Tumor cells at tumor invasive front/border':'Front tumor','Tumor stroma at tumor invasive front/border':'Front stroma'};a['plot_group']=a.Layer3.map(mapping).fillna('Other tissue');a.loc[a.Layer3.isin(['Tumor cells in mucosa and submucosa','Tumor cells in muscularis propria','Tumor cells in subserosa']),'plot_group']='Other tumor';a.loc[a.Layer3.isin(['Tumor stroma in mucosa and submucosa','Tumor stroma in muscularis propria','Tumor stroma in subserosa']),'plot_group']='Other tumor stroma';fig,ax=plt.subplots(figsize=(7,6),constrained_layout=True)
for lab,col,size in [('Other tissue','#DDDDDD',3),('Other tumor','#4F83B9',4),('Other tumor stroma','#7AAB87',4),('Front stroma','#A97242',9),('Front tumor','#D0473D',12)]:
 z=a[a.plot_group==lab];ax.scatter(z.pxl_col_in_fullres,z.pxl_row_in_fullres,c=col,s=size,label=lab,rasterized=True)
ax.invert_yaxis();ax.set_aspect('equal');ax.set_axis_off();ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False,markerscale=2);ax.set_title('M-ST-29 | Patient 11\nAuthor pathology labels on spatial spot coordinates',fontsize=12);fig.savefig(F/'route2_pathology_map.png',dpi=180);fig.savefig(F/'route2_pathology_map.svg');plt.close(fig)
(O/'results/INTEGRATION_COMPLETE.json').write_text(json.dumps({'exit_code':0,'figures':3,'route2_screen_genes':len(ss),'source_matched_genes':len(src),'scope':'Cross-dataset localization and descriptive covariance only; no inferred ligand receptor causality or validated novel targets.'},indent=2));print(out.sort_values('abs_effect',ascending=False)[['gene','median_delta','positive','negative','median_fibroblast_minus_tumor','fibroblast_higher']].head(20).to_string(index=False));print(corr.to_string())
