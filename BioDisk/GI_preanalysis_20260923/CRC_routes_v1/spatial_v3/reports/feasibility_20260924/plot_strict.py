import pathlib,sys
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3');V=O/'background_validation/run_20260924T052010Z/V01';a=pd.read_csv(O/'results/screen/patient_effects.tsv',sep='\t');a=a[(a.scope=='baseline')&(a.gene=='CGNL1')&(a.method=='pseudobulk')&(a.contrast=='tumor_near300_vs_far600')];b=pd.read_csv(V/'patient_effects.tsv',sep='\t');b=b[(b.gene=='CGNL1')&(b.method=='pseudobulk')&(b.rule=='5_of_6')];m=a.merge(b,on='patient',suffixes=('_initial','_strict'));out=O/'reports/feasibility_20260924';out.mkdir(parents=True,exist_ok=True)
fig,ax=plt.subplots(figsize=(8,4.2))
for _,r in m.iterrows():ax.plot([0,1],[r.effect_initial,r.effect_strict],'-o',color='#79919E',alpha=.8,markersize=4)
ax.scatter([0,1],[m.effect_initial.mean(),m.effect_strict.mean()],s=100,c='#C86B45',zorder=5,marker='D',label='Mean across the same 10 patients');ax.axhline(0,c='#344E5C',lw=.8);ax.set_xticks([0,1],['Initial spatial contrast','Mixed-edge exclusion']);ax.set_ylabel('Patient difference in log2(CPM+1)');ax.set_title('CGNL1: comparison restricted to the same 10 patients');ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,loc='upper right');fig.tight_layout();fig.savefig(out/'strict_matched_patients.png',dpi=180)
