from common import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Q/'CLOSEOUT'
p=pd.read_csv(Q/'F02_SPATIAL/patient_effects.tsv',sep='\t')
q=pd.read_csv(Q/'F04_MARKER_SENS/donor_state_effects.tsv',sep='\t')
# Final state and cross-cohort effect graphics; no inferential error bars for n3.
fig,axs=plt.subplots(1,2,figsize=(11.5,4.6))
src=pd.read_csv(Q/'F01_SOURCE_STATE/donor_state_effects.tsv',sep='\t');src=src[(src.lineage=='endothelial')&(src.target=='CGNL1')]
panels=['actomyosin','endothelial_junction','focal_adhesion','matrix_remodeling'];labels=['Actomyosin','EC junction','Focal adhesion','Matrix remodeling','EC junction\n(no label markers)']
for j,pan in enumerate(panels+['endothelial_junction_no_label_markers']):
 z=(q if j==4 else src[src.panel==pan]).effect.to_numpy();axs[0].scatter(np.arange(len(z))*.07-.07+j,z,s=36,color='#137F85');axs[0].plot([j-.16,j+.16],[np.median(z)]*2,color='#143A4B',lw=2)
axs[0].axhline(0,color='gray',lw=.7);axs[0].set_xticks(range(5),labels,rotation=20,ha='right');axs[0].set_ylabel('Matched mean log1p(CP10K) difference');axs[0].set_title('CGNL1 endothelial state | 3 donor labels')
main=p[np.isclose(p.q,.25)&np.isclose(p.pseudocount,1)]
for j,(g,cohort) in enumerate([('CGNL1','discovery'),('CGNL1','external_author_processed'),('PCSK5','discovery'),('PCSK5','external_author_processed')]):
 z=main[(main.gene==g)&(main.dataset==cohort)].effect.to_numpy();offset=np.linspace(-.13,.13,len(z));axs[1].scatter(j+offset,z,s=30,color='#137F85' if g=='CGNL1' else '#C4724E');axs[1].plot([j-.18,j+.18],[np.median(z)]*2,color='black',lw=2)
axs[1].axhline(0,color='gray',lw=.7);axs[1].set_xticks(range(4),['CGNL1\ndiscovery n10','CGNL1\nexternal n3','PCSK5\ndiscovery n10','PCSK5\nexternal n3']);axs[1].set_ylabel('Spatial difference in log2(CPM+1)');axs[1].set_title('Distance quartiles | same definition')
for ax in axs:ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(D/'final_evidence_comparison.png',dpi=190);plt.close(fig)
print('Final numerical audit and figure complete')
