import pathlib,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2/candidate_routes');F=O/'figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
genes=['TMEM91','CDH23','COPZ2'];colors=['#3a6e8b','#bc6b40','#6d8558','#876e9c']
p=pd.read_csv(O/'within_type_paired_effects.tsv',sep='\t');p=p[(p.level=='broad_type')&(p.min_cells==20)&(p.celltype=='Macrophage')];fig,ax=plt.subplots(1,3,figsize=(10,3.3),sharey=True)
for a,g in zip(ax,genes):
 z=p[p.gene==g].sort_values('patient')
 for i,(_,r) in enumerate(z.iterrows()):a.plot([0,1],[r.primary_log2CPM,r.ascites_log2CPM],'-o',c=colors[i],lw=1.4,label=r.patient.replace('MDA_',''))
 a.set_xticks([0,1],['Primary','Ascites']);a.set_title(g);a.set_xlim(-.15,1.15);a.grid(axis='y',alpha=.2)
ax[0].set_ylabel('Macrophage pseudobulk log2(CPM + 1)');ax[-1].legend(loc='best',fontsize=8,frameon=False);fig.suptitle('Paired macrophage expression: 4 independent patients',fontsize=12);fig.tight_layout();fig.savefig(F/'F1_paired_expression.png');fig.savefig(F/'F1_paired_expression.pdf');plt.close(fig)
p=pd.read_csv(O/'macrophage_subtype_fractions.tsv',sep='\t');fig,ax=plt.subplots(1,4,figsize=(10,3.3),sharey=True)
for a,ct in zip(ax,sorted(p.subtype.unique())):
 z=p[p.subtype==ct].sort_values('patient')
 for i,(_,r) in enumerate(z.iterrows()):a.plot([0,1],np.array([r.primary_fraction_in_macrophages,r.ascites_fraction_in_macrophages])*100,'-o',c=colors[i],lw=1.4)
 a.set_xticks([0,1],['Primary','Ascites']);a.set_title(ct);a.set_xlim(-.15,1.15);a.set_ylim(-3,103);a.grid(axis='y',alpha=.2)
ax[0].set_ylabel('Share of captured macrophages (%)');fig.suptitle('Author-defined macrophage subtype replacement: 4 paired patients',fontsize=12);fig.tight_layout();fig.savefig(F/'F2_subtype_composition.png');fig.savefig(F/'F2_subtype_composition.pdf');plt.close(fig)
x=pd.read_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t');scores=pd.read_csv(O/'external_program_scores.tsv',sep='\t');stats=pd.read_csv(O/'external_state_program_correlations.tsv',sep='\t');fig,ax=plt.subplots(1,3,figsize=(11,3.65),sharey=True)
for a,g in zip(ax,genes):
 z=x[(x.celltype=='Macrophage_enriched')&(x.gene==g)&(x.n_cells>=20)].merge(scores[(scores.celltype=='Macrophage_enriched')&(scores.program=='phagolysosomal')][['patient','score']],on='patient');st=stats[(stats.celltype=='Macrophage_enriched')&(stats.gene==g)&(stats.program=='phagolysosomal')].iloc[0]
 for typ,c in [('Malignant Ascites Fluid','#bc6b40'),('Peritoneal Washing Fluid','#3a6e8b')]:
  zz=z[z.sample_type==typ]
  if len(zz):a.scatter(zz.log2CPM,zz.score,c=c,label=typ.replace(' Fluid',''),s=28,alpha=.85)
 # Preserve any differently named type.
 for typ in set(z.sample_type)-{'Malignant Ascites Fluid','Peritoneal Washing Fluid'}:
  zz=z[z.sample_type==typ];a.scatter(zz.log2CPM,zz.score,label=typ,s=28)
 a.set_title(g);a.set_ylim(7.6,10.3);a.set_xlabel('Pseudobulk log2(CPM + 1)');a.text(.03,.97,f"n={len(z)}, rho={st.rho:.3f}\nP={st.P:.2g}, FDR={st.FDR_54:.3g}\nAdjusted rho={st.partial_rho:.3f}\nadjusted FDR={st.FDR_partial_54:.3g}",transform=a.transAxes,va='top',fontsize=10,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'});a.grid(alpha=.15)
ax[0].set_ylabel('Phagolysosomal marker score');ax[-1].legend(loc='lower right',fontsize=9,frameon=False);fig.suptitle('Independent cohort: patient-level state associations within macrophage-enriched cells',fontsize=11);fig.tight_layout();fig.savefig(F/'F3_external_state_correlations.png');fig.savefig(F/'F3_external_state_correlations.pdf');plt.close(fig)
print('Created 3 figure sets')

