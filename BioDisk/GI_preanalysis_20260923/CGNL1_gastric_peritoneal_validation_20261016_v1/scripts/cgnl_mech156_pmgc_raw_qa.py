from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,zipfile,gzip,io,csv,datetime,pandas as pd,numpy as np
O=Path.cwd()/'mechanism_investigation_v2';panels=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text())['source_panels'];spec=json.loads((O/'PMGC2026_QC_SPECIFICATION_145.json').read_text());wanted=set(sum(panels.values(),[])+['CGNL1']);hb=set(spec['primary_scRNA_QC']['erythrocyte_gene_set']);profiles=pd.read_csv(O/'results/PMGC2026_scRNA_source_profiles.tsv',sep='\t');records=[]
with zipfile.ZipFile(O/'raw/PMGC2026/Human_scRNA.zip') as z:
 for sample in ['PM1','PM2','PM3','PM4']:
  prefix='Human_scRNA/'+sample+'/'
  def rd(name):return io.TextIOWrapper(gzip.GzipFile(fileobj=z.open(prefix+name)),encoding='utf-8')
  with rd('features.tsv.gz') as f:ft=list(csv.reader(f,delimiter='\t'))
  with rd('barcodes.tsv.gz') as f:bar=[r[0] for r in csv.reader(f,delimiter='\t')]
  ix={i+1:r[1] for i,r in enumerate(ft) if r[1] in wanted};mt={i+1 for i,r in enumerate(ft) if r[1].startswith('MT-')};hi={i+1 for i,r in enumerate(ft) if r[1] in hb};N=len(bar);tot=np.zeros(N,dtype=np.int64);ng=np.zeros(N,dtype=np.int64);mito=np.zeros(N,dtype=np.int64);ery=np.zeros(N,dtype=np.int64);counts={g:np.zeros(N,dtype=np.int64) for g in wanted};nnz=0
  with rd('matrix.mtx.gz') as f:
   for line in f:
    if not line.startswith('%'):shape=list(map(int,line.split()));break
   assert shape[:2]==[len(ft),N]
   for line in f:
    a,b,c=line.split();i=int(a);j=int(b)-1;v=int(c);assert v>0;tot[j]+=v;ng[j]+=1;nnz+=1
    if i in mt:mito[j]+=v
    if i in hi:ery[j]+=v
    if i in ix:counts[ix[i]][j]+=v
  assert nnz==shape[2];d=pd.read_csv(O/'processed/PMGC2026_scRNA'/(sample+'_cell_source_QC.tsv.gz'),sep='\t').set_index('barcode').loc[bar];assert np.array_equal(d.total_UMI,tot) and np.array_equal(d.detected_genes,ng)
  for g in wanted:assert np.array_equal(d[g],counts[g]),(sample,g)
  mf=mito/np.maximum(tot,1);ef=ery/np.maximum(tot,1);assert np.max(abs(mf-d.mitochondrial_fraction))<1e-12 and np.max(abs(ef-d.erythrocyte_fraction))<1e-12;c=spec['primary_scRNA_QC'];primary=(ng>=c['genes_min_inclusive'])&(ng<=c['genes_max_inclusive'])&(tot<c['UMI_max_exclusive'])&(mf<c['mitochondrial_fraction_max_exclusive'])&(ef<c['erythrocyte_fraction_max_exclusive']);sensitivity=(tot>=500)&(ng>=200)&(ng<=7000)&(mf<=.25);assert np.array_equal(primary,d.primary_QC) and np.array_equal(sensitivity,d.sensitivity_QC)
  immune=sum((counts[g]>0).astype(int) for g in panels['immune'])>=2;meso=((counts['UPK3B']+counts['LRRN4'])>0)&((counts['WT1']+counts['CALB2']+counts['MSLN'])>0)&~immune;assert np.array_equal(meso,d.group.eq('strict_mesothelial'))
  for scheme,keep in [('primary',primary),('sensitivity',sensitivity)]:
   for row in profiles[(profiles['sample']==sample)&(profiles.scheme==scheme)].itertuples():
    mask=keep&d.group.eq(row.group).to_numpy();assert int(mask.sum())==row.n_cells and int(tot[mask].sum())==row.total_UMI and int(counts['CGNL1'][mask].sum())==row.CGNL1_UMI;assert int((counts['CGNL1'][mask]>0).sum())==row.CGNL1_positive_cells
    if row.total_UMI:assert abs(np.log1p(row.CGNL1_UMI/row.total_UMI*1e4)-row.CGNL1_ln1pCP10k)<1e-12
  records.append({'sample':sample,'triplets':nnz,'cells':N,'integer_target_marker_and_library_counts_equal':True,'QC_equal':True,'mesothelial_mask_equal':True,'source_aggregation_equal':True});print('RAW_QA',sample,nnz,flush=True)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'method':'independent line-by-line MatrixMarket counts versus scipy-based percell source export; QC and source totals recomputed','samples':records,'identity_limit':'numeric validity does not establish cell lineage or malignancy'};(O/'results/PMGC2026_RAW_NUMERIC_QA_156.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
