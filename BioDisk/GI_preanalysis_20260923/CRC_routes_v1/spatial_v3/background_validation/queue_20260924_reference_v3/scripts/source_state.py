from common import *
from scipy.stats import rankdata
from scipy.spatial import cKDTree
D=Q/'V06_SOURCE_STATE';D.mkdir(exist_ok=True);old=pathlib.Path(C['previous_queue'])/'V02_SOURCE';samples=json.loads((I/'SAMPLES.json').read_text());meta=pd.read_csv(old/'profile_metadata.tsv',sep='\t')
selected=set(meta[(meta.label=='epithelial_CNV_supported')&(meta.cells>=50)&meta.role.eq('liver_tumor')]['sample'])
panels={'junction':['TJP1','TJP2','OCLN','CDH1','CTNNA1','CTNNB1','CLDN3','CLDN4','CLDN7'],'actomyosin':['MYH9','MYH10','MYL9','MYLK','ROCK1','ROCK2'],'hepatic':['ALB','APOA1','APOA2','TTR','CPS1','HNF4A'],'cycling':['MKI67','TOP2A','CDK1','PCNA','MCM2','MCM5']}
dump(D/'PLAN.json',{'panels':panels,'unit':'donor labels','primary':'within-CNV-supported cells CGNL1-positive vs negative matched 1:1 without replacement on log1p UMI within0.25 and log1p detected genes within0.25; >=10 matched pairs/sample','endpoint':'mean difference log1p(CP10K) panel mean; target excluded; donor-equal sample averaging','scope':'Detection-associated state, not induction or causation. n<10 no inferential tests. Zero-positive samples retained in eligibility.'})
effects=[];qc=[]
for sample in samples:
 name=sample['sample']
 if name not in selected:continue
 X,sy,bc,lib=load_sample(sample);cl=pd.read_csv(old/'cells'/(name+'.tsv.gz'),sep='\t').set_index('barcode').reindex(bc);assert cl.label.notna().all()
 epi=(cl.label=='epithelial_CNV_supported').to_numpy();ix=np.flatnonzero(sy=='CGNL1');assert len(ix)==1;cg=X[ix[0]].toarray().ravel();ng=X.getnnz(axis=0);pos=np.flatnonzero(epi&(cg>0));neg=np.flatnonzero(epi&(cg==0))
 q={'sample':name,'dataset':sample['dataset'],'donor':sample['patient_label_provisional'],'CNV_cells':int(epi.sum()),'CGNL1_positive':len(pos),'positive_pct':float(len(pos)/max(epi.sum(),1)*100)}
 features=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for ip in sorted(pos,key=lambda i:(cg[i],i)):
  cand=np.array([j for j in neg if j not in used],dtype=int)
  if len(cand)==0:break
  dis=np.abs(features[cand]-features[ip]);ok=(dis<=.25).all(axis=1)
  if not ok.any():continue
  cand=cand[ok];delta=dis[ok];j=cand[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((ip,j))
 q['matched_pairs']=len(pairs);q['eligible']=len(pairs)>=10;qc.append(q)
 if len(pairs)>=10:
  p,n=np.array(pairs).T
  for panel,genes in panels.items():
   ids=np.flatnonzero(np.isin(sy,genes))
   if len(ids)<len(genes)*.7:continue
   A=np.log1p(X[ids][:,np.r_[p,n]].toarray()/lib[np.r_[p,n]]*1e4).mean(axis=0);diff=A[:len(p)]-A[len(p):];effects.append({**q,'panel':panel,'effect':float(diff.mean()),'genes_measured':len(ids),'mean_UMI_log_difference':float((features[p,0]-features[n,0]).mean()),'mean_gene_log_difference':float((features[p,1]-features[n,1]).mean())})
 pd.DataFrame(qc).to_csv(D/'eligibility.tsv',sep='\t',index=False)
pd.DataFrame(effects).to_csv(D/'sample_state_effects.tsv',sep='\t',index=False)
rows=[]
if effects:
 e=pd.DataFrame(effects);d=e.groupby(['dataset','donor','panel'],as_index=False).effect.mean();d.to_csv(D/'donor_state_effects.tsv',sep='\t',index=False)
 for panel,v in d.groupby('panel'):rows.append({'panel':panel,'donors':len(v),'positive':int((v.effect>0).sum()),'mean':float(v.effect.mean()),'median':float(v.effect.median())})
dump(D/'RESULT.json',{'status':'exploratory_state_computed' if effects else 'insufficient_positive_cells_for_matched_state','eligibility':qc,'summary':rows,'no_inferential_tests':True})
