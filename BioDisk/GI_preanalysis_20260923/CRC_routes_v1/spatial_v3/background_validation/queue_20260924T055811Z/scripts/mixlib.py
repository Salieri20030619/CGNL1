from common import *
from scipy.optimize import nnls
from scipy import stats

def fit_reference(X,sy,lib):
 r=np.load(Q/'V02_SOURCE/reference.npz');genes=r['genes'];types=r['types'];P=r['CPM'];unique=~pd.Index(sy).duplicated(keep=False);avail=set(sy[unique]);allowed=np.array([g in avail and g not in C['targets'] and not g.startswith(('MT-','RPL','RPS')) for g in genes]);chosen=set()
 for k in range(len(types)):
  other=(P.sum(axis=1)-P[:,k])/max(len(types)-1,1);rank=np.log2((P[:,k]+1)/(other+1));rank[~allowed|(P[:,k]<5)]=-np.inf;ii=np.argsort(rank)[::-1];chosen.update(ii[np.isfinite(rank[ii])][:C['mixture']['marker_top_per_type']].tolist())
 gi=np.array(sorted(chosen));assert len(gi)>=C['mixture']['min_fitting_genes'];index={g:i for i,g in enumerate(sy) if unique[i]};ix=np.array([index[g] for g in genes[gi]]);assert (ix>=0).all();weight=1/np.sqrt(P[gi].mean(axis=1)+1);A=P[gi]*weight[:,None];condition=float(np.linalg.cond(A));assert condition<C['mixture']['max_profile_condition'],'reference profile collinearity';coeff=np.zeros((X.shape[1],len(types)));validation=[]
 for start in range(0,X.shape[1],500):
  stop=min(start+500,X.shape[1]);Y=X[ix,start:stop].toarray()/np.maximum(lib[start:stop],1)*1e6;Y*=weight[:,None]
  for k in range(Y.shape[1]):coeff[start+k]=nnls(A,Y[:,k],maxiter=10*len(types))[0]
 frac=coeff/np.maximum(coeff.sum(axis=1,keepdims=True),1e-10);pred={}
 for gene in C['targets']:
  jj=np.flatnonzero(genes==gene)
  if len(jj)==1:pred[gene]=coeff@P[jj[0]]
 return types,frac,pred,{'fit_genes':len(gi),'profile_condition':condition,'targets_excluded_from_fitting':C['targets'],'measurement':'relative RNA mixture weights, not validated cell fractions'}
def summarize(e,path,family_n=4):
 if not e:dump(path/'RESULT.json',{'status':'insufficient_data'});return
 e=pd.DataFrame(e);e.to_csv(path/'section_effects.tsv',sep='\t',index=False);p=e.groupby(['patient','gene','mask']).effect.mean().dropna().reset_index();p.to_csv(path/'patient_effects.tsv',sep='\t',index=False);rows=[]
 for key,v in p.groupby(['gene','mask']):
  x=v.effect.values;n=len(x);r={'gene':key[0],'mask':key[1],'n_patients':n,'positive':int((x>0).sum()),'mean_effect':float(x.mean()),'P':np.nan,'CI_low':np.nan,'CI_high':np.nan}
  if n>=10:
   se=stats.sem(x);t=stats.t.ppf(.975,n-1);r.update(P=float(stats.ttest_1samp(x,0).pvalue),CI_low=float(x.mean()-t*se),CI_high=float(x.mean()+t*se))
  rows.append(r)
 s=pd.DataFrame(rows);s['FDR']=np.nan
 for mask,z in s.groupby('mask'):
  q=bh(z.P.values,4);q[z.P.isna().values]=np.nan;s.loc[z.index,'FDR']=q
 s.to_csv(path/'summary.tsv',sep='\t',index=False);dump(path/'RESULT.json',{'status':'exploratory_computed_pending_review','effects':rows,'interpretation':'Predicted-composition residuals do not establish malignant-cell causation.'})
