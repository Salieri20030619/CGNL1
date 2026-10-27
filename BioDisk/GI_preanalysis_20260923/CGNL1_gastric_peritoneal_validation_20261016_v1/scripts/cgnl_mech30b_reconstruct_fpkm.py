from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,hashlib,numpy as np,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2'
amend={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Author per-sample FPKM gzip fails decompression; annotated read counts and group FPKM pass full gzip CRC','method':'Reconstruct FPKM = read count / sum of gene counts per sample * 1e9 / author gene_length; validate all group means against intact author group-FPKM matrix before analysis','no_silent_replacement':True,'prior_outcome_knowledge':'Target CPM sensitivity already evaluated; this is normalization reconstruction validation, not a newly independent test','acceptance':'all matched genes/groups equal within absolute1e-6 or relative1e-7; otherwise do not label author-scale reconstructed'}
p=O/'FPKM_RECONSTRUCTION_AMENDMENT.json';assert p.exists()
r=pd.read_csv(O/'raw/GSE223165_readcount_genename.txt.gz',sep='\t',index_col=0,low_memory=False);m=pd.read_csv(O/'results/GSE223165_CPM_target_samples.tsv',sep='\t',dtype={'donor':str});cs=m.matrix_column.tolist();x=r[cs].astype(float);length=pd.to_numeric(r.gene_length);assert (length>0).all();z=x.div(x.sum(0),axis=1).div(length,axis=0)*1e9
g=pd.read_csv(O/'raw/GSE223165_fpkm_group.txt.gz',sep='\t',index_col=0);assert g.index.is_unique and z.index.is_unique and set(g.index)==set(z.index);z=z.loc[g.index];groupmap={'Control':'control','CM':'CM','ASC1':'ASC1','ASC2':'ASC2','ASC3':'ASC3'};audit=[]
for a,b in groupmap.items():
 rebuilt=z[m.loc[m.condition==b,'matrix_column']].mean(1);delta=abs(rebuilt-g[a]);okay=np.isclose(rebuilt,g[a],rtol=1e-7,atol=1e-6);audit.append({'group':a,'rows':len(g),'pass_rows':int(okay.sum()),'max_absolute_difference':float(delta.max()),'mean_absolute_difference':float(delta.mean())})
passed=all(a['pass_rows']==a['rows'] for a in audit);out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':passed,'audit':audit,'unit':'FPKM reconstructed from author counts and gene lengths, validated against author group values'};(O/'results/FPKM_RECONSTRUCTION_QA.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
if passed:z.to_csv(O/'raw/reconstructed_GSE223165_fpkm_sample.txt.gz',sep='\t',encoding='utf-8')
