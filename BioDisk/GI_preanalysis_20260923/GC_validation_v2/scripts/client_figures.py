import pathlib,json,numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');O=R/'04_results';F=O/'figures';F.mkdir(exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc');plt.rcParams.update({'font.family':'Noto Sans CJK JP','font.size':14,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
blue='#246B91';orange='#CB7040';grey='#7C8792'
def save(name):
 plt.tight_layout(pad=1.5);plt.savefig(F/(name+'.png'),dpi=220,bbox_inches='tight');plt.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close()
# Figure 1: all important negative sensitivity results shown.
z=pd.read_csv(O/'bulk_robustness'/'composition_proxy_sensitivity.tsv',sep='\t');loo=pd.read_csv(O/'bulk_robustness'/'leave_one_gene_out.tsv',sep='\t');main=pd.read_csv(O/'bulk_robustness'/'primary_effect_CI.tsv',sep='\t').set_index('metric').value
labels=['主模型','去除 FABP4','调整脂肪细胞标记','调整脂肪与髓系标记'];vals=[main['effect.Estimate'],loo.loc[loo.omitted_gene=='FABP4','effect.Estimate'].iloc[0],z.loc[z.model=='adipocyte','effect.Estimate'].iloc[0],z.loc[z.model=='adipocyte_myeloid','effect.Estimate'].iloc[0]];los=[main['lo.Estimate'],loo.loc[loo.omitted_gene=='FABP4','lo.Estimate'].iloc[0],z.loc[z.model=='adipocyte','lo.Estimate'].iloc[0],z.loc[z.model=='adipocyte_myeloid','lo.Estimate'].iloc[0]];his=[main['hi.Estimate'],loo.loc[loo.omitted_gene=='FABP4','hi.Estimate'].iloc[0],z.loc[z.model=='adipocyte','hi.Estimate'].iloc[0],z.loc[z.model=='adipocyte_myeloid','hi.Estimate'].iloc[0]]
fig,ax=plt.subplots(1,2,figsize=(11,4));y=np.arange(4);ax[0].errorbar(vals,y,xerr=[np.array(vals)-los,np.array(his)-vals],fmt='o',color=blue,capsize=3);ax[0].axvline(0,c=grey,lw=1);ax[0].set_yticks(y,labels);ax[0].invert_yaxis();ax[0].set_xlabel('腹膜转移 − 原发灶\n脂质面板差值及95%置信区间');ax[0].set_title('A  组成因素解释了较多组织差异',loc='left',fontsize=12)
d=pd.read_csv(O/'bulk_robustness'/'paired_program_differences.tsv',sep='\t');ax[1].scatter(d.adipocyte,d.lipid,s=40,c=blue);ax[1].set_xlabel('同患者脂肪细胞标记变化');ax[1].set_ylabel('同患者脂质面板变化');ax[1].set_title('B  变化方向与富脂环境相连',loc='left',fontsize=12);ax[1].text(.03,.97,'17对患者\nSpearman ρ = 0.931\nFDR = 2.26 × 10^-7',va='top',transform=ax[1].transAxes,fontsize=10);save('F1_tissue_source')
# Figure 2: paired within macrophages, no inferential test in four patients.
p=pd.read_csv(O/'paired_scRNA'/'paired_gene_effects.tsv',sep='\t');p=p[p.celltype=='Macrophage'];fig,ax=plt.subplots(1,2,figsize=(9,3.8))
for a,g in zip(ax,['CD36','TREM2']):
 q=p[p.gene==g].sort_values('patient')
 for i,row in enumerate(q.itertuples()):a.plot([0,1],[row.primary_log2CPM,row.ascites_log2CPM],'o-',lw=1.3,alpha=.8,label=row.patient.replace('MDA_',''))
 a.set_xticks([0,1],['原发灶巨噬细胞','腹水巨噬细胞']);a.set_ylabel('患者聚合 log2(CPM + 1)');a.set_title(f'{g}  4/4对升高',fontsize=13);a.set_xlim(-.2,1.2)
ax[1].legend(fontsize=8,frameon=False);save('F2_paired_macrophages')
# Figure 3: independent localization; no metastatic vs nonmetastatic claim.
x=pd.read_csv(O/'independent_scRNA_refined'/'marker_expression_patient_compartment.tsv',sep='\t');x=x[x.n_cells>=20];fig,ax=plt.subplots(1,2,figsize=(9,4))
for a,g in zip(ax,['CD36','TREM2']):
 w=x[(x.gene==g)&x.celltype.isin(['Macrophage_enriched','Lymphoid'])].pivot(index='patient',columns='celltype',values='log2CPM').dropna()
 for row in w.itertuples():a.plot([0,1],[row.Lymphoid,row.Macrophage_enriched],'-',c='#BBCBD6',lw=.7,alpha=.7)
 a.scatter(np.zeros(len(w)),w.Lymphoid,c=grey,s=13);a.scatter(np.ones(len(w)),w.Macrophage_enriched,c=blue,s=13);a.set_xticks([0,1],['淋巴细胞群','C1Q巨噬细胞富集群']);a.set_ylabel('患者聚合 log2(CPM + 1)');a.set_title(f'{g}  独立28例一致定位',fontsize=13);a.set_xlim(-.2,1.2);a.text(.03,.97,'配对P = 7.45 × 10^-9\n3终点FDR = 1.12 × 10^-8',va='top',transform=a.transAxes,fontsize=9)
save('F3_independent_localization')
# Figure 4 for internal appendix: whole marker source, low support types explicitly omitted.
c=pd.read_csv(O/'paired_scRNA'/'cell_source_patient_summary.tsv',sep='\t');c=c[(c.patients>=3)&c.gene.isin(['CD36','TREM2','FABP4','LPL','APOE','GPNMB','SPP1'])];w=c.pivot(index='celltype',columns='gene',values='median_patient_log2CPM');fig,ax=plt.subplots(figsize=(8,3));im=ax.imshow(w.values,cmap='Blues');ax.set_yticks(range(len(w)),w.index);ax.set_xticks(range(w.shape[1]),w.columns);fig.colorbar(im,ax=ax,label='患者中位 log2(CPM + 1)');ax.set_title('参考单细胞来源 仅显示至少3位患者覆盖的类型');save('S1_cell_source')
print(json.dumps({'figures':sorted(p.name for p in F.glob('*.png'))},ensure_ascii=False))
