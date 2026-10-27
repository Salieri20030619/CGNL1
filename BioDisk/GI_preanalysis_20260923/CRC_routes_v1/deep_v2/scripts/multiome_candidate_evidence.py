import pathlib,json,gzip,re,numpy as np,pandas as pd,h5py
from scipy import sparse
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');O=B/'CRC_unused'
V=B/'CRC_routes_v1/deep_v2';R=V/'results/multiome';R.mkdir(exist_ok=True);genes=json.loads((V/'PLAN.json').read_text())['frozen_genes'];coords={}
with gzip.open(V/'metadata/GRCh38.p14.ncbi.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if c[2]!='gene':continue
  a=re.fullmatch(r'NC_0000(\d{2})\.\d+',c[0])
  if not a or int(a.group(1))>24:continue
  g=re.search(r'(?:^|; )gene "([^"]+)"',c[8]);g=g.group(1) if g else re.search(r'gene_id "([^"]+)"',c[8]).group(1)
  if g in [x.replace('.','-') for x in genes]:
   chrom='chr'+({23:'X',24:'Y'}.get(int(a.group(1)),str(int(a.group(1)))));tss=int(c[3])-1 if c[6]=='+' else int(c[4])-1;coords.setdefault(g,[]).append((chrom,max(0,tss-(2000 if c[6]=='+' else 500)),tss+(500 if c[6]=='+' else 2000),c[6]))
rows=[];peaks=[];qc=[]
for name,kind in [('GSE294559_TP001.h5','tumor_tissue'),('GSE294559_TO001.h5','tumor_organoid')]:
 with h5py.File(B/'02_raw'/name,'r') as f:
  z=f['matrix'];X=sparse.csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=tuple(z['shape'][:]));sy=np.array([x.decode() for x in z['features']['name'][:]]);ty=np.array([x.decode() for x in z['features']['feature_type'][:]]);ri=np.flatnonzero(ty=='Gene Expression');pi=np.flatnonzero(ty=='Peaks');RNA=X[ri].tocsr();AT=X[pi].tocsr();rs=sy[ri];ps=sy[pi];rnaUMI=np.asarray(RNA.sum(axis=0)).ravel();atUMI=np.asarray(AT.sum(axis=0)).ravel();rg=np.asarray((RNA>0).sum(axis=0)).ravel();ap=np.asarray((AT>0).sum(axis=0)).ravel();mt=np.asarray(RNA[np.char.startswith(rs.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(rnaUMI,1);keep=(rnaUMI>=500)&(rg>=200)&(ap>=200)&(mt<=.25);ep=np.asarray((RNA[np.isin(rs,['EPCAM','KRT8','KRT18','KRT19'])]>0).sum(axis=0)).ravel()>=2;epi=keep&ep
  # Parse called ATAC peak intervals; compare promoter evidence qualitatively, not cross-assay count fold changes.
  pc=np.array([x.split(':')[0] for x in ps]);starts=np.array([int(x.split(':')[1].split('-')[0]) for x in ps]);ends=np.array([int(x.split(':')[1].split('-')[1]) for x in ps]);qc.append({'sample':name,'kind':kind,'patient':'P001','barcodes':X.shape[1],'jointQC':int(keep.sum()),'epithelial_like_jointQC':int(epi.sum())})
  for gene in genes:
   g=gene.replace('.','-');gi=np.flatnonzero(rs==g);co=coords.get(g,[])
   if len(gi)!=1 or len(co)!=1:continue
   chrom,st,en,strand=co[0];ix=np.flatnonzero((pc==chrom)&(starts<en)&(ends>st));rv=np.asarray(RNA[gi].sum(axis=0)).ravel();av=np.asarray(AT[ix].sum(axis=0)).ravel() if len(ix) else np.zeros(X.shape[1]);peaks.append({'gene':gene,'sample':name,'chromosome':chrom,'window_start':st,'window_end':en,'called_peaks':';'.join(ps[ix])})
   for pop,mask in [('all_jointQC',keep),('epithelial_like_jointQC',epi)]:
    n=int(mask.sum());rows.append({'gene':gene,'patient':'P001','kind':kind,'population':pop,'n_cells':n,'promoter_overlapping_called_peaks':len(ix),'RNA_detected_cells':int((rv[mask]>0).sum()),'ATAC_promoter_detected_cells':int((av[mask]>0).sum()),'both_detected_cells':int(((rv[mask]>0)&(av[mask]>0)).sum()),'RNA_detection_fraction':float((rv[mask]>0).mean()),'promoter_detection_fraction':float((av[mask]>0).mean()),'RNA_pseudobulk_CPM':float(rv[mask].sum()/max(rnaUMI[mask].sum(),1)*1e6),'ATAC_promoter_per10k_peakcounts':float(av[mask].sum()/max(atUMI[mask].sum(),1)*1e4)})
a=pd.DataFrame(rows);a.to_csv(R/'candidate_RNA_ATAC_evidence.tsv',sep='\t',index=False);pd.DataFrame(peaks).to_csv(R/'candidate_promoter_windows.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(R/'QC.tsv',sep='\t',index=False);(R/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'independent_patients':1,'samples':2,'scope':'Same-cell RNA/ATAC promoter-window feasibility only; no PT-LM contrast; no inferential tests; called peak sets differ; gene-boundary TSS proxy not enhancer causality'},indent=2));print(pd.DataFrame(qc).to_string(index=False));print(a[(a.gene.isin(['CGNL1','CREB3L1','PCSK5','SLC9B2']))&(a.population=='epithelial_like_jointQC')].to_string(index=False))
