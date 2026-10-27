from common import *
from scipy.stats import spearmanr
D=Q/'V02_SOURCE';D.mkdir(exist_ok=True);old=pathlib.Path(C['previous_queue'])/'V02_SOURCE';m=pd.read_csv(old/'profile_metadata.tsv',sep='\t');keep=m[(m.cells>=50)&m.role.str.startswith('liver')]
zhep=np.load(Q/'NORMAL_REF/hepatocyte_profiles.npz');cache={};allgenes=set(zhep['genes'])
for name in keep['sample'].unique():
 z=np.load(old/'profiles'/(name+'.npz'));cache[name]=z;allgenes&=set(z['genes'])
genes=np.array(sorted(allgenes));assert len(genes)>10000
types=C['reference_gate']['required']+['epithelial_unresolved','biliary_epithelial_like'];bytype={};source=[]
for ct in types:
 vals={}
 for donor,v in keep[keep.label==ct].groupby('donor_label'):
  a=[]
  for _,r in v.iterrows():
   z=cache[r['sample']];a.append(z['profiles'][pd.Index(z['genes']).get_indexer(genes),np.flatnonzero(z['labels']==ct)[0]])
  vals[donor]=np.mean(a,axis=0);source.append({'type':ct,'donor':donor,'origin':'previous_raw_count_reference'})
 if ct=='hepatocyte':
  gi=pd.Index(zhep['genes']).get_indexer(genes)
  for j,donor in enumerate(zhep['donors']):vals[donor]=zhep['CPM'][gi,j];source.append({'type':ct,'donor':donor,'origin':'normal_liver_processed_reference'})
 if len(vals)>=3:bytype[ct]=vals
required=C['reference_gate']['required'];gate={ct:len(bytype.get(ct,{})) for ct in required};allowed=all(v>=3 for v in gate.values());dump(D/'GATE.json',{'allowed':allowed,'donor_labels_per_required_type':gate,'reason':'unchanged >=3 donor and >=50 cell coverage; normal liver supplement','scope':'Exploratory cross-platform reference. Existing cancer-cohort donor labels provisional.'});assert allowed
types=list(bytype);P=np.stack([np.mean(list(bytype[ct].values()),axis=0) for ct in types],axis=1)
np.savez_compressed(D/'reference.npz',genes=genes,types=np.array(types),CPM=P)
Ps=P.copy();Ps[:,types.index('hepatocyte')]=np.mean([v for k,v in bytype['hepatocyte'].items() if k.startswith('GSE115469:')],axis=0)
np.savez_compressed(D/'reference_normal_hep_only.npz',genes=genes,types=np.array(types),CPM=Ps)
pd.DataFrame(source).to_csv(D/'reference_donors.tsv',sep='\t',index=False)
audit=[]
for donor,v in bytype['hepatocyte'].items():
 oth=np.mean([x for k,x in bytype['hepatocyte'].items() if k!=donor],axis=0);select=~np.isin(genes,C['targets'])&~np.char.startswith(genes,'MT-')
 audit.append({'donor':donor,'Spearman_to_other_hepatocyte_donors':float(spearmanr(np.log1p(v[select]),np.log1p(oth[select])).statistic),**{g:float(v[np.flatnonzero(genes==g)[0]]) for g in C['targets'] if g in genes}})
pd.DataFrame(audit).to_csv(D/'hepatocyte_reference_agreement.tsv',sep='\t',index=False)
dump(D/'RESULT.json',{'status':'reference_ready_pending_fit_QC','common_genes':len(genes),'types':types,'gate':gate,'sensitivity':'replace merged hepatic profile with external normal liver only; other types unchanged','normalization':C['normal_reference_scale']})
