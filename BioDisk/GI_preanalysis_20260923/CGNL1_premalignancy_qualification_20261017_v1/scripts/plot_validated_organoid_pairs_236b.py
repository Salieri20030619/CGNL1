
from pathlib import Path
import json
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path.cwd();out=p/'results/GSE210993_bulk_230';figdir=out/'figures';figdir.mkdir(exist_ok=True)
assert json.loads((out/'INDEPENDENT_NUMERICAL_QA_232.json').read_text())['pass']
meta=pd.DataFrame(json.loads((p/'metadata/GSE210993_ANALYSIS_SAMPLES_227.json').read_text()))
target=pd.read_csv(out/'target_samples.tsv',sep=chr(9))
patients=json.loads((p/'STATISTICAL_SPECIFICATION_227.json').read_text())['paired_patient_ids']
keep=target[~target.reclassified&target.patient_id.isin(patients)&target.type.isin(['Human IM organoid','Human normal gastric organoid'])]
w=keep.groupby(['patient_id','type']).CGNL1_log2TPM1.mean().unstack().loc[patients]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={'width_ratios':[1,1.45]})
colors=plt.cm.tab10(np.linspace(0,0.8,len(patients)))
for patient,color in zip(patients,colors):
 y=[w.loc[patient,'Human normal gastric organoid'],w.loc[patient,'Human IM organoid']]
 axes[0].plot([0,1],y,'o-',color=color,linewidth=1.2,markersize=5,label=patient)
axes[0].set(xticks=[0,1],xticklabels=['Normal gastric','IM'],ylabel='Patient mean log2(TPM + 1)',title='CGNL1: 7 paired patients')
axes[0].text(.02,.03,'7/7 lower in IM; exact paired P = 0.015625',transform=axes[0].transAxes,fontsize=9)
axes[0].legend(frameon=False,fontsize=8,loc='upper right',ncol=2)
context=pd.read_csv(out/'context_gene_samples.tsv',sep=chr(9)).merge(meta[['sample_name','patient_id','type','reclassified']],on='sample_name',validate='many_to_one')
genes=['CGNL1','CGN','CAMSAP3','CDH1','TJP1','EPCAM']
source=[]
for k,g in enumerate(genes):
 if g=='CGNL1':d=(w['Human IM organoid']-w['Human normal gastric organoid']).to_numpy()
 else:
  a=context[(context.gene==g)&~context.reclassified&context.patient_id.isin(patients)&context.type.isin(['Human IM organoid','Human normal gastric organoid'])]
  a=a.groupby(['patient_id','type']).log2TPM1.mean().unstack().loc[patients]
  d=(a['Human IM organoid']-a['Human normal gastric organoid']).to_numpy()
 axes[1].scatter(k+np.linspace(-.15,.15,len(d)),d,s=24,color='#b23639' if g=='CGNL1' else '#326b95',alpha=.85)
 axes[1].plot([k-.23,k+.23],[np.median(d)]*2,color='black',linewidth=2)
 source.extend({'gene':g,'patient':pat,'difference':float(v)} for pat,v in zip(patients,d))
axes[1].axhline(0,color='#888888',linewidth=.8)
axes[1].set(xticks=range(len(genes)),xticklabels=genes,ylabel='IM minus normal: delta log2(TPM + 1)',title='Predefined junction context: paired differences')
axes[1].tick_params(axis='x',labelrotation=30)
fig.text(.53,.015,'Each dot = one patient; bars = medians. Context genes descriptive.',fontsize=8)
for ax in axes:
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#eeeeee',linewidth=.6);ax.set_axisbelow(True)
fig.tight_layout(pad=1.8,rect=(0,.06,1,1));fig.savefig(figdir/'CGNL1_paired_organoid_and_junction_context_v2.png',dpi=200);fig.savefig(figdir/'CGNL1_paired_organoid_and_junction_context_v2.svg')
assert (figdir/'plot_source_values.tsv').exists()
print(str(figdir/'CGNL1_paired_organoid_and_junction_context_v2.png'))
