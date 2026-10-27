import pathlib,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2/pair_validation');F=O/'figures';F.mkdir(exist_ok=True);a=pd.read_csv(O/'tissue_gene_compartment_values.tsv.gz',sep='\t');qc=pd.read_csv(O/'tissue_sample_qc.tsv',sep='\t');raw=pd.read_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t');ps=pd.read_csv(O/'external_program_scores.tsv',sep='\t');tests=pd.read_csv(O/'external_state_program_correlations.tsv',sep='\t')
# Algebraic decomposition on captured cells only, sample-weighted equally within groups.
rows=[];summ=[]
for g in ['RAP1GAP2','OPLAH']:
 z=a[(a.gene==g)&~a.celltype.isin(['All_QC','Epithelial_strict'])].merge(qc[['sample','qc_cells']],on='sample');z['fraction']=z.n_cells/z.qc_cells;z['contribution']=z.fraction*z.mean_CP10K
 assert np.allclose(z.groupby('sample').fraction.sum(),1)
 total=a[(a.gene==g)&(a.celltype=='All_QC')].set_index('sample');assert np.allclose(z.groupby('sample').contribution.sum().sort_index(),total.mean_CP10K.sort_index())
 for ct in z.celltype.unique():
  h={}
  for group in ['PT','PM']:
   t=z[(z.celltype==ct)&(z.group==group)];p=t.fraction.sum()/3;c=t.contribution.sum()/3;h[group]=(p,c,c/p if p>0 else np.nan)
  p0,c0,m0=h['PT'];p1,c1,m1=h['PM'];state=(p0+p1)/2*(m1-m0);comp=(m0+m1)/2*(p1-p0);rows.append({'gene':g,'celltype':ct,'PT_mean_fraction':p0,'PM_mean_fraction':p1,'PT_conditional_mean_CP10K':m0,'PM_conditional_mean_CP10K':m1,'expression_term':state,'composition_term':comp,'total_contribution_change':c1-c0,'unresolved_if_absent':c1-c0-state-comp})
 rr=pd.DataFrame(rows);rr=rr[rr.gene==g];delta=total[total.group=='PM'].mean_CP10K.mean()-total[total.group=='PT'].mean_CP10K.mean();assert np.isclose(rr.total_contribution_change.sum(),delta);summ.append({'gene':g,'total_mean_CP10K_change':delta,'expression_terms_sum':rr.expression_term.sum(),'composition_terms_sum':rr.composition_term.sum(),'scope':'captured-cell algebra only; expression term includes subtype mixture and does not establish cell-intrinsic regulation'})
pd.DataFrame(rows).to_csv(O/'captured_expression_decomposition.tsv',sep='\t',index=False);pd.DataFrame(summ).to_csv(O/'captured_expression_decomposition_summary.tsv',sep='\t',index=False)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(9,7.1),constrained_layout=True)
for i,g in enumerate(['RAP1GAP2','OPLAH']):
 for j,ct in enumerate(['Epithelial_like','Epithelial_strict']):
  ax=axs[i,j];z=a[(a.gene==g)&(a.celltype==ct)]
  for k,(group,color) in enumerate([('PT','#28649A'),('PM','#CA663E')]):
   t=z[z.group==group].sort_values('sample');off=np.linspace(-.09,.09,len(t));ax.scatter(k+off,t.log2CPM,c=color,s=55,zorder=3)
   for x,(_,v) in zip(k+off,t.iterrows()):
    left=v['sample'].endswith(('Ca_2','F_14'));ax.annotate(v['sample'].split('_',1)[1],(x,v.log2CPM),xytext=(-6 if left else 6,5),ha='right' if left else 'left',textcoords='offset points',fontsize=8)
   ax.hlines(t.log2CPM.median(),k-.22,k+.22,color=color,lw=2)
  ax.set_xticks([0,1],['Primary (3 samples)','Peritoneal (3 samples)']);ax.set_xlim(-.4,1.65);ax.set_title(g+' | '+('Broad epithelial-like' if j==0 else 'Strict epithelial gate'));ax.set_ylabel('Pseudobulk log2(CPM + 1)');ax.grid(axis='y',alpha=.15)
fig.suptitle('Independent tissue validation: GSE308231\nEach point is one specimen; descriptive only, no paired donor mapping',fontsize=12)
fig.savefig(F/'tissue_candidate_validation.png',dpi=200);fig.savefig(F/'tissue_candidate_validation.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(9,4.4),constrained_layout=True)
for ax,(g,program) in zip(axs,[('RAP1GAP2','adhesion'),('OPLAH','glutathione_cycle')]):
 z=raw[(raw.gene==g)&(raw.celltype=='Epithelial_like')&(raw.n_cells>=20)].merge(ps[(ps.celltype=='Epithelial_like')&(ps.program==program)][['patient','score']],on='patient');r=tests[(tests.gene==g)&(tests.celltype=='Epithelial_like')&(tests.program==program)].iloc[0]
 for (t,v),col in zip(z.groupby('sample_type'),['#28649A','#CA663E']):ax.scatter(v.log2CPM,v.score,s=32,label=t,color=col,alpha=.85)
 ax.set_xlabel(g+' log2(CPM + 1)');ax.set_ylabel(program.replace('_',' ')+' score');ax.set_title(f'{g}: n={len(z)} patients\nRaw rho={r.rho:.3f}, FDR={r.FDR_24:.3g}\nAdjusted rho={r.partial_rho:.3f}, FDR={r.FDR_partial_24:.3g}',fontsize=10);ax.legend(fontsize=8,loc='best',frameon=False);ax.grid(alpha=.15)
fig.suptitle('GSE228598 | independent patient-level associations\nBH within 24 tests per model; adjustment: sample type, cell count, median UMI',fontsize=11);fig.savefig(F/'external_program_associations.png',dpi=200);fig.savefig(F/'external_program_associations.svg');plt.close(fig)
# Audit shared state marker profile and strict-gate representation without relabeling based on candidate expression.
C=O/'epithelial_states';c=pd.read_csv(C/'cell_clusters.tsv.gz',sep='\t');strict=c[c.strict_epithelial==True].groupby(['sample','group','cluster']).size().rename('strict_cells').reset_index();strict.to_csv(C/'strict_gate_cluster_counts.tsv',sep='\t',index=False)
m=pd.read_csv(C/'state_sample_expression.tsv.gz',sep='\t');m[(m.cluster==1)&(m.n_cells>=20)&m.gene.isin(['EPCAM','KRT8','KRT18','KRT19','PTPRC','LST1','COL1A1','VIM','MUC1','MSLN','CALB2'])].groupby(['group','gene']).agg(n_samples=('sample','nunique'),median_log2CPM=('log2CPM','median'),median_detection=('detection_fraction','median')).reset_index().to_csv(C/'shared_cluster1_marker_audit.tsv',sep='\t',index=False)
(O/'FIGURE_DECOMPOSITION_COMPLETE.json').write_text(json.dumps({'exit_code':0,'figures':2,'decomposition_identity_checked':True,'scope':'descriptive, algebraic captured-cell decomposition not causal attribution'},indent=2));print(pd.DataFrame(summ).to_string(index=False));print(pd.DataFrame(rows).query("gene=='OPLAH' and celltype=='Epithelial_like'").to_string(index=False));print(strict[strict.cluster==1].to_string(index=False))
