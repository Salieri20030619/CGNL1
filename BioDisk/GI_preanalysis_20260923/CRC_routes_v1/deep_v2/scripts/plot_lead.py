import pathlib,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');F=O/'results/figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
m=pd.read_csv(O/'results/bulk_validation/sample_values.tsv',sep='\t');e=pd.read_csv(O/'results/bulk_validation/log2FPKM_filtered.tsv.gz',sep='\t',index_col=0);z=m.copy();z['expression']=e.loc['CGNL1',z['sample']].values;z=z.pivot(index='patient',columns='site',values='expression');a=pd.read_csv(O/'results/second_validation/focused_sample_values.tsv',sep='\t');a=a[a.gene=='CGNL1'];b=a[a.site.isin(['PT','LM'])&a.strict_colon].pivot(index='patient',columns='site',values='expression');fig,axs=plt.subplots(1,2,figsize=(9,4.2),constrained_layout=True)
for ax,data,title,ylab in [(axs[0],z,'GSE50760: 18 paired CRC patients','log2(FPKM + 1)'),(axs[1],b,'GSE14297: 13 strict-colon pairs','Author-normalized array expression')]:
 for _,r in data.iterrows():ax.plot([0,1],[r.PT,r.LM],color='#C6CBD1',lw=.8,zorder=1)
 ax.scatter(np.zeros(len(data)),data.PT,color='#2B6899',s=24,zorder=2);ax.scatter(np.ones(len(data)),data.LM,color='#C55B3B',s=24,zorder=2);ax.set_xticks([0,1],['Primary','Liver metastasis']);ax.set_xlim(-.2,1.2);ax.set_ylabel(ylab);ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.15)
axs[0].text(.03,.98,'16/18 increased\nP = 1.48e-4; genome FDR = 0.00532',transform=axs[0].transAxes,va='top',fontsize=9);axs[1].text(.03,.98,'12/13 increased\nP = 3.12e-6; genome FDR = 0.000997',transform=axs[1].transAxes,va='top',fontsize=9)
for ax in axs:lo,hi=ax.get_ylim();ax.set_ylim(lo,hi+(hi-lo)*.25)
fig.savefig(F/'lead_independent_replication.png',dpi=180);fig.savefig(F/'lead_independent_replication.svg');plt.close(fig)
c=pd.read_csv(O/'results/cnv/candidate_patient_effects.tsv',sep='\t');c=c[c.candidate=='CGNL1'];fig,axs=plt.subplots(1,3,figsize=(10.5,3.8),constrained_layout=True)
for ax,sub,title in [(axs[0],'CNV95','CNV support: 95th reference cutoff'),(axs[1],'CNV99','CNV support: 99th reference cutoff')]:
 t=c[c.subset==sub]
 for _,r in t.iterrows():ax.plot([0,1],[r.primary_log2CPM,r.liver_log2CPM],color='#BCC3CC',lw=1)
 ax.scatter(np.zeros(len(t)),t.primary_log2CPM,color='#2B6899',s=26);ax.scatter(np.ones(len(t)),t.liver_log2CPM,color='#C55B3B',s=26);ax.set_xticks([0,1],['Primary','Liver']);ax.set_ylabel('log2(CPM + 1)');ax.set_title(title+'\n4 paired patients, descriptive',fontsize=9)
ax=axs[2]
for i,(site,label,col) in enumerate([('PT','Primary','#2B6899'),('LM','Liver met.','#C55B3B'),('NL','Normal liver','#666666')]):
 t=a[(a.site==site)&(a.strict_colon if site!='NL' else True)].expression;ax.scatter(i+np.linspace(-.09,.09,len(t)),t,color=col,s=23);ax.hlines(t.median(),i-.18,i+.18,color=col,lw=2)
ax.set_xticks([0,1,2],['Primary\nn=13','Liver met.\nn=13','Normal liver\nn=5']);ax.set_title('Normal liver has higher expression\nSeparate controls; descriptive',fontsize=9);ax.set_ylabel('Array-normalized expression');fig.savefig(F/'lead_source_and_boundary.png',dpi=180);fig.savefig(F/'lead_source_and_boundary.svg');plt.close(fig)
profiles=json.loads((O/'results/cnv/group_profiles.json').read_text());z=np.array([r['profile'] for r in profiles]);labs=[r['patient']+' '+r['site']+' '+r['broad'] for r in profiles];fig,ax=plt.subplots(figsize=(11,8),constrained_layout=True);im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-.3,vmax=.3,interpolation='nearest');ax.set_yticks(range(len(labs)),labs,fontsize=7);ax.set_xlabel('Genomically ordered smoothed windows (autosomes)');ax.set_title('RNA-inferred CNV: patient / tissue / author lineage averages\nCandidates and hepatic marker genes excluded; not DNA confirmation',fontsize=11);fig.colorbar(im,ax=ax,label='Mean inferred deviation',shrink=.5);fig.savefig(F/'cnv_reference_audit.png',dpi=160);plt.close(fig)
(F/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'figures':3,'gene':'CGNL1'},indent=2))
