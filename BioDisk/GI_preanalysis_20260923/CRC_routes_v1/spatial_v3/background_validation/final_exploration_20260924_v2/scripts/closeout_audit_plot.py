from common import *
from scipy import stats
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Q/'CLOSEOUT';D.mkdir(exist_ok=True)
s=pd.read_csv(Q/'F02_SPATIAL/summary.tsv',sep='\t');p=pd.read_csv(Q/'F02_SPATIAL/patient_effects.tsv',sep='\t')
for _,r in s.iterrows():
 v=p[(p.dataset==r.dataset)&(p.gene==r.gene)&np.isclose(p.q,r.q)&np.isclose(p.pseudocount,r.pseudocount)].effect.to_numpy()
 if len(v)>=10:
  pv=stats.ttest_1samp(v,0).pvalue;se=stats.sem(v);w=stats.t.ppf(.975,len(v)-1)*se
  assert np.isclose(pv,r.P) and np.allclose([v.mean()-w,v.mean()+w],[r.CI_low,r.CI_high])
 else:assert pd.isna(r.P) and pd.isna(r.FDR)
for key,v in s.groupby(['dataset','q','pseudocount']):
 raw=np.nan_to_num(v.P.to_numpy(),nan=1);order=np.argsort(raw);adj=np.minimum.accumulate((raw[order]*2/np.arange(1,len(raw)+1))[::-1])[::-1];q=np.empty(len(raw));q[order]=np.minimum(adj,1);good=v.P.notna().to_numpy();assert np.allclose(q[good],v.FDR.to_numpy()[good])
for stage in ['F01_SOURCE_STATE','F04_MARKER_SENS']:
 e=pd.read_csv(Q/stage/'sample_state_effects.tsv',sep='\t');v=pd.read_csv(Q/stage/'donor_state_effects.tsv',sep='\t');re=e.groupby(['dataset','donor','lineage','target','panel'],as_index=False).effect.mean();m=v.merge(re,on=['dataset','donor','lineage','target','panel'],suffixes=('_saved','_re'));assert len(m)==len(v) and np.allclose(m.effect_saved,m.effect_re)
 assert (e.matched_pairs>=10).all() and (e.abs_logUMI_match_difference<=.25+1e-8).all() and (e.abs_logGene_match_difference<=.25+1e-8).all()
q=pd.read_csv(Q/'F04_MARKER_SENS/donor_state_effects.tsv',sep='\t');assert len(q)==3 and (q.effect>0).all()
exitstatus=json.loads((Q/'QUEUE_EXIT.json').read_text());assert exitstatus['exit_code']==0
dump(D/'AUDIT.json',{'status':'passed','spatial_CI_P_two_gene_FDR':'passed','source_donor_equal_aggregation':'passed','matching_caliper':'passed','external_small_n_no_inference':'passed','nonoverlap_endothelial_panel':'3/3 positive','primary_queue_exit_code':0,'marker_sensitivity_native_exit_code':0,'exploration_closed':True})
# Final state and cross-cohort effect graphics; no inferential error bars for n3.
fig,axs=plt.subplots(1,2,figsize=(11.5,4.6))
src=pd.read_csv(Q/'F01_SOURCE_STATE/donor_state_effects.tsv',sep='\t');src=src[(src.lineage=='endothelial')&(src.target=='CGNL1')]
panels=['actomyosin','endothelial_junction','focal_adhesion','matrix_remodeling'];labels=['Actomyosin','EC junction','Focal adhesion','Matrix remodeling','EC junction\n(no label markers)']
for j,pan in enumerate(panels+['endothelial_junction_no_label_markers']):
 z=(q if j==4 else src[src.panel==pan]).effect.to_numpy();axs[0].scatter(np.arange(len(z))*.07-.07+j,z,s=36,color='#137F85');axs[0].plot([j-.16,j+.16],[np.median(z)]*2,color='#143A4B',lw=2)
axs[0].axhline(0,color='gray',lw=.7);axs[0].set_xticks(range(5),labels,rotation=20,ha='right');axs[0].set_ylabel('Matched difference in mean log1p(CP10K)');axs[0].set_title('CGNL1 endothelial state | 3 donor labels')
main=p[np.isclose(p.q,.25)&np.isclose(p.pseudocount,1)]
for j,(g,cohort) in enumerate([('CGNL1','discovery'),('CGNL1','external_author_processed'),('PCSK5','discovery'),('PCSK5','external_author_processed')]):
 z=main[(main.gene==g)&(main.dataset==cohort)].effect.to_numpy();offset=np.linspace(-.13,.13,len(z));axs[1].scatter(j+offset,z,s=30,color='#137F85' if g=='CGNL1' else '#C4724E');axs[1].plot([j-.18,j+.18],[np.median(z)]*2,color='black',lw=2)
axs[1].axhline(0,color='gray',lw=.7);axs[1].set_xticks(range(4),['CGNL1\ndiscovery n10','CGNL1\nexternal n3','PCSK5\ndiscovery n10','PCSK5\nexternal n3']);axs[1].set_ylabel('Within-tumor distance-ranked log2(CPM+1) difference');axs[1].set_title('Same distance-quartile definition | descriptive comparison')
for ax in axs:ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(D/'final_evidence_comparison.png',dpi=190);plt.close(fig)
print('Final numerical audit and figure complete')
