from common import *
from scipy.stats import spearmanr
A=Q/'AUTHOR_VALIDATION';X=sparse.load_npz(A/'author_counts.npz');sy=pd.read_csv(A/'author_genes.tsv',sep='\t').gene.astype(str).to_numpy();v=pd.read_csv(A/'author_spots.tsv',sep='\t');bar=pd.read_csv(A/'author_barcodes.tsv',sep='\t');v['barcode']=bar.barcode;out=[]
for k in [4,5]:
 t=v[v.batch==k];Y=X[:,t.index];name='patient'+str(k);z=np.load(Q/'V04_INPUT_QC'/(name+'_features.npz'));Z=sparse.load_npz(Q/'V04_INPUT_QC'/(name+'_counts.npz'));bi=pd.Index(z['barcodes']).get_indexer(t.barcode);rg=pd.Index(z['genes']);unique=~rg.duplicated(keep=False);idx={g:j for j,g in enumerate(rg) if unique[j]};ids=np.array([i for i,g in enumerate(sy) if g in idx]);gi=np.array([idx[sy[i]] for i in ids]);B=Y[ids].tocsr();C0=Z[gi][:,bi].tocsr();d=B-C0;d.eliminate_zeros();rr,cc=d.nonzero();orig=np.asarray(C0[rr,cc]).ravel();new=np.asarray(B[rr,cc]).ravel();rec={'patient':k,'shape':B.shape,'mismatch_entries':d.nnz,'author_nonzero':B.nnz,'raw_nonzero':C0.nnz,'sum_author':float(B.sum()),'sum_raw':float(C0.sum()),'column_sum_correlation':float(spearmanr(np.asarray(B.sum(axis=0)).ravel(),np.asarray(C0.sum(axis=0)).ravel()).statistic),'examples':[{'gene':sy[ids[r]],'barcode':t.barcode.iloc[c],'author':float(n),'raw':float(o)} for r,c,n,o in zip(rr[:10],cc[:10],new[:10],orig[:10])]}
 for g in ['CGNL1','ALB','EPCAM']:
  a=Y[np.flatnonzero(sy==g)[0]].toarray().ravel();b=Z[idx[g],bi].toarray().ravel();rec[g]={'author_sum':float(a.sum()),'raw_sum':float(b.sum()),'equal_spots':int((a==b).sum()),'rho':float(spearmanr(a,b).statistic)}
 out.append(rec)
dump(A/'COUNT_DIAGNOSIS.json',out)
