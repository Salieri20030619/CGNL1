
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,gzip,csv,itertools,datetime,hashlib
import numpy as np
p=Path.cwd();m=p/'mechanism_investigation_v2'
spec=json.loads((m/'GSE286057_COMPARATOR_SPECIFICATION_195.json').read_text(encoding='utf-8'))
needed=sorted(set(sum(spec['source_panels'].values(),[])+['CGNL1','CGN','CAMSAP3','FN1','TNC','ACTA2','TAGLN','KRT8','KRT18','KRT19','PAX8','ITLN1']))
checks=[];out=m/'results/GSE286057_comparator_198'
assert json.loads((out/'STATISTICS_COMPLETE_204.json').read_text())['original_matrix_independent_QA_pending']
for s in spec['samples']:
 prefix=s['accession']+'_'+s['title'];root=m/'raw/GSE286057_comparator_195'
 with gzip.open(root/(prefix+'_features.tsv.gz'),'rt',encoding='utf-8') as f:features=list(csv.reader(f,delimiter=chr(9)))
 with gzip.open(root/(prefix+'_barcodes.tsv.gz'),'rt',encoding='utf-8') as f:barcodes=[r[0] for r in csv.reader(f,delimiter=chr(9))]
 with gzip.open(m/'processed/GSE286057_comparator_198'/(prefix+'_cells.tsv.gz'),'rt',encoding='utf-8') as f:cells=list(csv.DictReader(f,delimiter=chr(9)))
 assert barcodes==[r['barcode'] for r in cells] and len(set(barcodes))==len(barcodes)
 nc=len(cells);ng=len(features);slots={g:i for i,g in enumerate(needed)}
 gene_slot=np.array([slots.get(r[1],-1) for r in features],dtype=np.int64)
 mitochondrial=np.array([r[1].startswith('MT-') for r in features])
 umi=np.zeros(nc,dtype=np.int64);detected=np.zeros(nc,dtype=np.int64);mt=np.zeros(nc,dtype=np.int64);target=np.zeros((len(needed),nc),dtype=np.int64)
 seen=0
 with gzip.open(root/(prefix+'_matrix.mtx.gz'),'rt',encoding='utf-8') as f:
  for line in f:
   if not line.startswith('%') and line.strip():dims=list(map(int,line.split()));break
  assert dims[:2]==[ng,nc]
  while True:
   lines=list(itertools.islice(f,250000))
   if not lines:break
   lines=[line for line in lines if line.strip()]
   a=np.fromstring(''.join(lines),sep=' ',dtype=np.int64)
   assert a.size==3*len(lines)
   a=a.reshape(-1,3);i=a[:,0]-1;j=a[:,1]-1;v=a[:,2]
   assert (i>=0).all() and (i<ng).all() and (j>=0).all() and (j<nc).all() and (v>=0).all()
   umi+=np.bincount(j,weights=v,minlength=nc).astype(np.int64)
   detected+=np.bincount(j[v>0],minlength=nc).astype(np.int64)
   mask=mitochondrial[i];mt+=np.bincount(j[mask],weights=v[mask],minlength=nc).astype(np.int64)
   slot=gene_slot[i];mask=slot>=0;np.add.at(target,(slot[mask],j[mask]),v[mask]);seen+=len(lines)
 assert seen==dims[2]
 def original(key):return np.array([float(r[key]) for r in cells])
 assert np.array_equal(umi,original('UMI')) and np.array_equal(detected,original('n_genes'))
 error=float(np.max(np.abs(mt/np.maximum(umi,1)-original('mito_fraction'))));assert error<1e-12
 for k,g in enumerate(needed):assert np.array_equal(target[k],original(g)),(s['title'],g)
 for scheme in ['primary','sensitivity']:
  q=spec['qc_'+scheme];mask=(umi>=q['UMI_min'])&(detected>=q['genes_min'])&(detected<=q['genes_max'])&(mt/np.maximum(umi,1)<=q['mitochondrial_fraction_max'])
  stored=np.array([r['qc_'+scheme]=='True' for r in cells]);assert np.array_equal(mask,stored)
 checks.append({'sample':s['title'],'cells':nc,'matrix_triplets':seen,'genes_checked':len(needed),'all_integer_counts_exact':True,'mito_max_error':error})
 (out/'RAW_MATRIX_QA_PROGRESS_211.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
 print('RAW_QA_PASS',s['title'],seen,flush=True)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'samples':len(checks),'method':'Independent streaming MatrixMarket text parser (no scipy mmread), raw integer accumulation; barcodes, every cell library size/detected genes/mitochondrial fraction, all target and marker counts, both QC sets compared exactly','checks':checks,'does_not_establish':'host genotype, clinical causation or novel mechanism','existing_results_unchanged':True}
(out/'RAW_MATRIX_INDEPENDENT_QA_211.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('ALL_36_RAW_MATRICES_QA_PASS',flush=True)
