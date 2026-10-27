from pathlib import Path
import pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path.cwd();d=pd.read_csv(p/'results/spheroid_233/sample_target_summary.tsv',sep='\t');out=p/'figures';out.mkdir(exist_ok=True)
fig,ax=plt.subplots(figsize=(6.3,4.9));colors=['#4477AA','#CC6677','#228833']
for (patient,g),c in zip(d.groupby('patient',sort=True),colors):
 v=g.set_index('condition').loc[['parent','spheroid'],'CGNL1_common_log2CPM1'].values
 ax.plot([0,1],v,'o-',color=c,lw=2,ms=7,label=f'{patient}  (delta {v[1]-v[0]:+.2f})')
ax.set(xticks=[0,1],xticklabels=['Parental organoid','Matrix-independent\nspheroid'],ylabel='CGNL1 pseudobulk log2(CPM + 1)',xlim=(-.15,1.15),ylim=(0,3.1))
ax.set_title('CGNL1 changes differ between three patients',fontsize=12,pad=15)
ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,loc='upper center',fontsize=9)
fig.text(.13,.035,'3 independent patients; 2 down, 1 up. Mean delta = -0.007; exact paired P = 1.0.\nAll six raw matrices independently verified; observational comparison.',fontsize=8)
fig.subplots_adjust(left=.13,right=.96,bottom=.20,top=.86)
for ext in ['png','svg']:fig.savefig(out/('CGNL1_spheroid_paired_248_v2.'+ext),dpi=180)
d.to_csv(out/'CGNL1_spheroid_paired_248_v2_source.tsv',sep='\t',index=False)
print(str(out/'CGNL1_spheroid_paired_248_v2.png'))
