from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,re,hashlib,datetime,itertools
import pandas as pd,numpy as np
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2';spec=O/'EXPOSURE_ANALYSIS_SPECIFICATION.json';assert spec.exists()
meta=(O/'metadata/GSE223165_samples.soft').read_bytes().decode('utf-8');rows=[]
for b in meta.split('^SAMPLE = ')[1:]:
 ls=b.splitlines()
 def field(k):return [l.split(' = ',1)[1] for l in ls if l.startswith('!Sample_'+k+' = ')]
 title=field('title')[0];donor=re.search(r'patient (\d+)',title).group(1);original=field('description')[0];col=original.replace(' ','_')
 condition='control' if 'untreated' in title else ('CM' if 'conditioned media' in title else 'ASC'+re.search(r'ascites sample(\d)',title).group(1))
 rows.append({'accession':ls[0],'title':title,'donor':donor,'description_original':original,'matrix_column':col,'condition':condition})
m=pd.DataFrame(rows);assert len(m)==41 and m.donor.nunique()==10 and not m[['donor','condition']].duplicated().any()
assert json.loads((O/'results/FPKM_RECONSTRUCTION_QA.json').read_text())['all_pass'];p=O/'raw/reconstructed_GSE223165_fpkm_sample.txt.gz';x=pd.read_csv(p,sep='\t',index_col=0);assert list(x.columns)==m.matrix_column.tolist();assert x.index.is_unique;assert np.isfinite(x.values).all() and (x.values>=0).all()
target='ENSG00000128849';assert target in x.index
a=pd.read_csv(O/'raw/GSE223165_readcount_genename.txt.gz',sep='\t',index_col=0,low_memory=False);assert a.loc[target,'gene_name']=='CGNL1'
m['FPKM']=x.loc[target,m.matrix_column].to_numpy();m['log2FPKM1']=np.log2(m.FPKM+1);m.to_csv(O/'results/GSE223165_reconstructedFPKM_target_samples.tsv',sep='\t',index=False,encoding='utf-8')
v=m.pivot(index='donor',columns='condition',values='log2FPKM1');delta=pd.DataFrame(index=v.index);delta['CM_minus_control']=v.CM-v.control;delta['mean_ASC_minus_control']=v[['ASC1','ASC2','ASC3']].mean(axis=1)-v.control
for c in ['ASC1','ASC2','ASC3']:delta[c+'_minus_control']=v[c]-v.control
delta.to_csv(O/'results/GSE223165_reconstructedFPKM_target_donor_deltas.tsv',sep='\t',encoding='utf-8')
def exact(d):
 d=d[d!=0];r=stats.rankdata(abs(d));obs=abs(np.sum(np.sign(d)*r));vals=[abs(np.dot(s,r)) for s in itertools.product([-1,1],repeat=len(d))];return float(np.mean(np.array(vals)>=obs-1e-12)) if len(d) else 1.
out=[]
for c in ['CM_minus_control','mean_ASC_minus_control']:
 d=delta[c].dropna().values;n=len(d);se=stats.sem(d);mu=np.mean(d);q=stats.t.ppf(.975,n-1)*se
 out.append({'contrast':c,'n_independent_mesothelial_donors':n,'mean_log2FPKM1_difference':mu,'median_difference':np.median(d),'positive':int((d>0).sum()),'negative':int((d<0).sum()),'zero':int((d==0).sum()),'primary_signedrank_P':exact(d),'secondary_paired_t_P':stats.ttest_1samp(d,0).pvalue,'mean_CI_low':mu-q,'mean_CI_high':mu+q,'min_leave_one_out_mean':min(np.delete(d,i).mean() for i in range(n)),'max_leave_one_out_mean':max(np.delete(d,i).mean() for i in range(n))})
r=pd.DataFrame(out);r['primary_FDR_two_target_contrasts']=stats.false_discovery_control(r.primary_signedrank_P);r.to_csv(O/'results/GSE223165_reconstructedFPKM_target_contrasts.tsv',sep='\t',index=False,encoding='utf-8')
print(r.to_string(index=False));print(delta.to_string());print(m.groupby('condition').FPKM.agg(['count','mean','median']).to_string())
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spec_sha256':hashlib.sha256(spec.read_bytes()).hexdigest(),'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'shape':list(x.shape),'sample_mapping':'GEO original description spaces explicitly mapped to matrix underscores; exact 41-column identity; donor/treatment from titles','source':'ovarian cancer conditioned media / ascites exposure to human primary mesothelial cells','mechanism_established':False,'gastric_specific_validation':False,'analysis_role':'prespecified FPKM scale reconstructed from intact counts and author lengths, all58735genes in5group means validated against author FPKM; not direct intact sample-FPKM download','counts_normalization_sensitivity':'results/GSE223165_CPM_target_contrasts.tsv'}
(O/'results/EXPOSURE_RECONSTRUCTED_FPKM_COMPLETE.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
