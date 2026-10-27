from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,hashlib,itertools
import pandas as pd,numpy as np
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2'
panels={'mesothelial_identity':['UPK3B','LRRN4','WT1','CALB2','MSLN'],'mesothelial_specific':['UPK3B','LRRN4'],'junction':['CDH1','TJP1','OCLN','JUP','CTNNB1'],'mesenchymal':['VIM','FN1','COL1A1','COL1A2','ACTA2','TAGLN'],'cycling':['MKI67','TOP2A','PCNA','UBE2C','CENPF']}
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Distinguish target decline from parallel lineage marker decline, not causal adjustment','posthoc':True,'known_prior_result':'CGNL1 declines after ovarian ascites; CM not nominal','panels':panels,'measure':'mean log2(FPKM+1) across fixed mapped unique genes; target-minus-identity score is relative-expression contrast, not regression-mediated effect','families':{'target_relative_identity':'2identity panels x2exposures=4','state_programs':'3programs x2exposures=6'},'independent_unit':'mesothelial donor; ASC3 averaged within donor','coverage':'all fixed genes retained when uniquely mapped; ambiguous/missing reported; zero-baseline identity panel invalidates interpretation; per-gene donor deltas descriptive','qualifiers':'n7 exploratory even if FDR passes; ovarian exposure not gastric, no protein activity inference'}
p=O/'IDENTITY_COMPARATOR_SPECIFICATION.json';assert not p.exists();p.write_text(json.dumps(spec,indent=2),encoding='utf-8')
x=pd.read_csv(O/'raw/reconstructed_GSE223165_fpkm_sample.txt.gz',sep='\t',index_col=0);a=pd.read_csv(O/'raw/GSE223165_readcount_genename.txt.gz',sep='\t',index_col=0,low_memory=False);m=pd.read_csv(O/'results/GSE223165_reconstructedFPKM_target_samples.tsv',sep='\t',dtype={'donor':str});y=np.log2(x+1);scores=m[['donor','condition','matrix_column']].copy();coverage=[];genes=set(sum(panels.values(),[])+['CGNL1','ANGPTL4','STC1']);gvalues={}
for g in sorted(genes):
 ids=a.index[a.gene_name==g].tolist()
 if len(ids)==1:gvalues[g]=y.loc[ids[0],scores.matrix_column].to_numpy()
 coverage.append({'gene':g,'n_rows':len(ids),'control_median_FPKM':float(x.loc[ids[0],m.loc[m.condition=='control','matrix_column']].median()) if len(ids)==1 else None})
for panel,gs in panels.items():
 good=[g for g in gs if g in gvalues];assert len(good)>=len(gs)*.8;scores[panel]=np.mean([gvalues[g] for g in good],axis=0)
scores['CGNL1']=gvalues['CGNL1']
for panel in ['mesothelial_identity','mesothelial_specific']:scores['CGNL1_relative_'+panel]=scores.CGNL1-scores[panel]
def exact(d):
 d=d[d!=0];r=stats.rankdata(abs(d));obs=abs(np.dot(np.sign(d),r));return float(np.mean([abs(np.dot(s,r))>=obs-1e-12 for s in itertools.product([-1,1],repeat=len(d))])) if len(d) else 1.
out=[];deltas=[]
for measure in ['CGNL1_relative_mesothelial_identity','CGNL1_relative_mesothelial_specific','junction','mesenchymal','cycling']:
 v=scores.pivot(index='donor',columns='condition',values=measure)
 for exposure in ['CM','mean_ASC']:
  d=(v.CM-v.control if exposure=='CM' else v[['ASC1','ASC2','ASC3']].mean(1)-v.control).dropna();out.append({'measure':measure,'exposure':exposure,'family':'target_relative_identity' if measure.startswith('CGNL1_relative') else 'state_programs','n_donors':len(d),'mean_difference':d.mean(),'median_difference':d.median(),'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P':exact(d.to_numpy())});deltas.extend({'measure':measure,'exposure':exposure,'donor':donor,'difference':val} for donor,val in d.items())
r=pd.DataFrame(out);r['FDR']=np.nan
for fam,ix in r.groupby('family').groups.items():r.loc[ix,'FDR']=stats.false_discovery_control(r.loc[ix,'P'])
gene_rows=[]
for g,vv in gvalues.items():
 t=scores[['donor','condition']].copy();t['value']=vv;v=t.pivot(index='donor',columns='condition',values='value')
 for exposure in ['CM','mean_ASC']:
  d=(v.CM-v.control if exposure=='CM' else v[['ASC1','ASC2','ASC3']].mean(1)-v.control).dropna();gene_rows.append({'gene':g,'exposure':exposure,'n_donors':len(d),'mean_difference':d.mean(),'positive':int((d>0).sum()),'negative':int((d<0).sum())})
for name,df in [('identity_comparator_tests',r),('identity_comparator_deltas',pd.DataFrame(deltas)),('identity_comparator_scores',scores),('identity_comparator_gene_coverage',pd.DataFrame(coverage)),('identity_comparator_gene_directions',pd.DataFrame(gene_rows))]:df.to_csv(O/'results'/('GSE223165_'+name+'.tsv'),sep='\t',index=False,encoding='utf-8')
print(r.to_string(index=False));print('COVERAGE',pd.DataFrame(coverage).to_string(index=False));print('IDENTITY_DIRECTIONS',pd.DataFrame(gene_rows).query("gene in ['CGNL1','UPK3B','LRRN4','WT1','CALB2','MSLN','ANGPTL4','STC1']").to_string(index=False))
