from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/annotation_panel_20261007_v1')
review=json.loads((R/'review/cluster_annotation_review.json').read_text());mapping={k:v['label'] for k,v in review['clusters'].items()}
def read(name):
 t=pd.read_csv(R/'results'/name,sep='\t',dtype=str);return t.set_index(t.columns[0]).astype(float)
m=read('refined_cluster_mean_log.tsv');f=read('refined_cluster_detection_fraction.tsv')
genes='HNF1A G6PC CYP2E1 APOB FGB FLT1 CLEC14A PTPRB ROBO4 PLVAP CD34 RGS5 MYH11 ACTA2 TAGLN COL1A1 COL1A2 DCN PDGFRA EPCAM KRT19 KRT7 CD3D CD3E GNLY NKG7 MS4A1 CD19 CD79A MZB1 TNFRSF17 CD68 CD163 C1QA MERTK LYZ CD1C CLEC9A TPSB2 KIT FCGR3B CSF3R CEACAM8'.split()
z=(m[genes]-m[genes].mean())/m[genes].std().replace(0,1)
fig,ax=plt.subplots(figsize=(19,11))
for j,g in enumerate(genes):artist=ax.scatter(np.repeat(j,len(m)),np.arange(len(m)),s=f[g].to_numpy()*95,c=z[g].clip(-2,2),cmap='RdBu_r',vmin=-2,vmax=2)
ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=9);ax.set_yticks(range(len(m)),[f'{cl}: {mapping[cl]}' for cl in m.index],fontsize=9);ax.invert_yaxis();ax.set_title('Panel marker evidence for reviewed expression clusters\nBefore per-cell marker-support and mixed-lineage filtering',fontsize=13,pad=14)
cb=fig.colorbar(artist,ax=ax,shrink=.32,pad=.02);cb.set_label('Mean log expression\nstandardized within each gene',fontsize=10)
handles=[ax.scatter([],[],s=v*95,color='#777777',label=f'{v:.0%}') for v in [.25,.5,1.]]
ax.legend(handles=handles,title='Fraction detected',loc='upper center',bbox_to_anchor=(.5,-.14),ncol=3,frameon=False)
fig.tight_layout();fig.savefig(R/'figures/marker_dotplot.png',dpi=180);plt.close(fig)
print('MARKER_REVIEW_PLOT_COMPLETE')
