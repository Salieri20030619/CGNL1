from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,tarfile,gzip,io,csv,numpy as np,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';spec=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text());panels=spec['source_panels'];known=pd.read_csv(O/'results/GSE308231_source_profiles.tsv',sep='\t');needed=set(sum(panels.values(),[])+['CGNL1']);intent={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'role':'posthoc identity robustness and independent raw triplet audit; aggregate CGNL1 outcomes already known','primary_source_results_unchanged':True,'additional_marker_rules':['specific_total_UMI>=2 and support_total_UMI>=2','bothUPK3B_and_LRRN4_detected and support>=1'],'not_identity_proof':'marker strength does not establish nonmalignancy; no cell-level inferential P','QC_sets_not_nested':'primary has no genes_max; sensitivity excludes>7000genes; inspect rather than assume additionalcells mustaddtargetUMI'};(O/'TISSUE_IDENTITY_AUDIT_SPEC_151.json').write_text(json.dumps(intent,indent=2),encoding='utf-8');cellrows=[];summary=[];qa=[]
with tarfile.open(O/'raw/GSE308231_RAW.tar') as tf:
 for sample in spec['samples']:
  acc=sample['accession'];prefix=next(n[:-len('_matrix.mtx.gz')] for n in tf.getnames() if n.startswith(acc) and n.endswith('_matrix.mtx.gz'))
  def text(suffix):return io.TextIOWrapper(gzip.GzipFile(fileobj=tf.extractfile(prefix+suffix)),encoding='utf-8')
  with text('_features.tsv.gz') as f:features=list(csv.reader(f,delimiter='\t'))
  with text('_barcodes.tsv.gz') as f:bar=[r[0] for r in csv.reader(f,delimiter='\t')]
  chosen={i+1:r[1] for i,r in enumerate(features) if r[1] in needed};assert len(chosen)==len(needed);assert [r[0] for r in features if r[1]=='CGNL1']==['ENSG00000128849'];N=len(bar);counts={g:np.zeros(N,dtype=np.int64) for g in needed};tot=np.zeros(N,dtype=np.int64);ng=np.zeros(N,dtype=np.int64);mt=np.zeros(N,dtype=np.int64);mi={i+1 for i,r in enumerate(features) if r[1].startswith('MT-')};nnz=0
  with text('_matrix.mtx.gz') as f:
   for line in f:
    if not line.startswith('%'):shape=list(map(int,line.split()));break
   assert shape[:2]==[len(features),N]
   for line in f:
    a,b,c=line.split();i=int(a);col=int(b)-1;v=int(c);assert v>0;tot[col]+=v;ng[col]+=1;nnz+=1
    if i in mi:mt[col]+=v
    if i in chosen:counts[chosen[i]][col]+=v
  assert nnz==shape[2];mito=mt/np.maximum(tot,1)
  def det(gs):return sum((counts[g]>0).astype(int) for g in gs)
  immune=det(panels['immune'])>=2;specific=sum(counts[g] for g in panels['mesothelial_specific']);support=sum(counts[g] for g in panels['mesothelial_support']);meso=(specific>0)&(support>0)&~immune;epi=(counts['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)&(specific==0)&~immune;schemes={}
  for scheme in ['primary','sensitivity']:
   q=spec['qc_'+scheme];keep=(tot>=q['UMI_min'])&(ng>=q['genes_min'])&(mito<=q['mitochondrial_fraction_max'])
   if 'genes_max' in q:keep&=ng<=q['genes_max']
   schemes[scheme]=keep
   for label,mask in [('strict_mesothelial',meso),('conservative_epithelial',epi)]:
    mask=mask&keep;k=known[(known['sample']==acc)&(known.scheme==scheme)&(known.group==label)].iloc[0];assert int(mask.sum())==k.n_cells and int(counts['CGNL1'][mask].sum())==k.target_UMI;assert int(tot[mask].sum())==k.total_UMI
   for label,mask in [('original',meso),('marker_UMI2',meso&(specific>=2)&(support>=2)),('both_specific',meso&(counts['UPK3B']>0)&(counts['LRRN4']>0))]:
    mask&=keep;n=int(mask.sum());total=int(tot[mask].sum());target=int(counts['CGNL1'][mask].sum());summary.append({'sample':acc,'scheme':scheme,'identity_rule':label,'n_cells':n,'CGNL1_UMI':target,'total_UMI':total,'CGNL1_detection':float((counts['CGNL1'][mask]>0).mean()) if n else None,'CGNL1_ln1pCP10k':float(np.log1p(target/total*1e4)) if total else None,'meets20cell_coverage':n>=20})
  for i in np.where(meso&(schemes['primary']|schemes['sensitivity']))[0]:cellrows.append({'sample':acc,'barcode':bar[i],'total_UMI':int(tot[i]),'detected_genes':int(ng[i]),'mitochondrial_fraction':float(mito[i]),'primary_QC':bool(schemes['primary'][i]),'sensitivity_QC':bool(schemes['sensitivity'][i]),**{g:int(counts[g][i]) for g in sorted(needed)}})
  qa.append({'sample':acc,'triplets':nnz,'raw_cells':N,'mesothelial_and_epithelial_counts_equal':True});print('AUDIT',acc,'triplets',nnz,flush=True)
pd.DataFrame(cellrows).to_csv(O/'results/GSE308231_mesothelial_cell_audit.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(summary).to_csv(O/'results/GSE308231_identity_sensitivity_151.tsv',sep='\t',index=False,encoding='utf-8');(O/'results/TISSUE_RAW_COUNT_QA_151.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'method':'line-by-lineMatrixMarket integer accumulation independent of scipy mmread','samples':qa,'posthoc_identity_sensitivity':True},indent=2),encoding='utf-8');print(pd.DataFrame(summary).to_string(index=False),flush=True)
