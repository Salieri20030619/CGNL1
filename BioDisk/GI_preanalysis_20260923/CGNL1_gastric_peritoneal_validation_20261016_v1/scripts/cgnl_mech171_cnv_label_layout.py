from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,pandas as pd,numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path.cwd()/'mechanism_investigation_v2';D=O/'results/CNV_identity/GSE308231_GSM9240675';order_groups=['immune_like','endothelial_like','fibroblast_like','mural_like','overlapping_stromal','unclassified','strict_mesothelial','conservative_epithelial']
for model in ['primary','immune_reference']:
 C=sparse.load_npz(D/(model+'_RNA_CNV.npz'));obs=pd.read_csv(D/(model+'_cell_descriptors.tsv'),sep='\t');pos=json.loads((D/(model+'_genomic_windows.json')).read_text(encoding='utf-8'))['chr_pos'];ix=np.concatenate([np.flatnonzero(obs.group.eq(g)) for g in order_groups]);fig,ax=plt.subplots(figsize=(15,10));im=ax.imshow(C[ix].toarray(),aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-.5,vmax=.5);ys=[];labels=[];start=0
 for g in order_groups:
  n=int(obs.group.eq(g).sum())
  if n:ys.append(start+(n-1)/2);labels.append(g+' (n='+str(n)+')');start+=n;ax.axhline(start-.5,color='#333333',lw=.4)
 ax.set_yticks(ys,labels,fontsize=8);positions=sorted(pos.items(),key=lambda x:x[1]);ticks=[];chrom=[]
 for i,(ch,start) in enumerate(positions):
  end=positions[i+1][1] if i+1<len(positions) else C.shape[1];ticks.append((start+end-1)/2);chrom.append(ch);ax.axvline(start-.5,color='#333333',lw=.3)
 ax.set_xticks(ticks,chrom,rotation=90,fontsize=8);ax.set_title('GSE308231_GSM9240675 | '+model+' | RNA-derived chromosomal patterns',fontsize=12);ax.set_xlabel('Autosomal expression windows; not measured DNA copy numbers');fig.colorbar(im,ax=ax,label='Centered smoothed relative RNA expression',fraction=.025,pad=.02);fig.tight_layout();fig.savefig(D/(model+'_heatmap_v2.png'),dpi=150);plt.close(fig)
print('Replotted2crowdedpanelsatgreaterheight; originalfigures andallnumericoutputs preserved')
