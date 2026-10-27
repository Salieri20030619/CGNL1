from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,tarfile,zipfile,gzip,hashlib,pandas as pd,numpy as np,anndata as ad
from scipy import io,sparse
from threadpoolctl import threadpool_limits
O=Path.cwd()/'mechanism_investigation_v2';spec=json.loads((O/'CNV_IDENTITY_SPECIFICATION_165.json').read_text(encoding='utf-8'));tissue=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text(encoding='utf-8'));panels=tissue['source_panels'];out=O/'processed/CNV_identity_inputs';out.mkdir(exist_ok=True);records=[]
for case_index,case in enumerate(spec['cases']):
 sample=case['sample'];dataset=case['dataset']
 if dataset=='GSE308231':
  archive=tarfile.open(O/'raw/GSE308231_RAW.tar');prefix=next(n[:-len('_matrix.mtx.gz')] for n in archive.getnames() if n.startswith(sample) and n.endswith('_matrix.mtx.gz'))
  def rd(name):return gzip.GzipFile(fileobj=archive.extractfile(prefix+'_'+name))
 else:
  archive=zipfile.ZipFile(O/'raw/PMGC2026/Human_scRNA.zip');prefix='Human_scRNA/'+sample+'/'
  def rd(name):return gzip.GzipFile(fileobj=archive.open(prefix+name))
 ft=pd.read_csv(rd('features.tsv.gz'),sep='\t',header=None);bar=pd.read_csv(rd('barcodes.tsv.gz'),header=None,sep='\t')[0].to_numpy(dtype=str);genes=ft[1].to_numpy(dtype=str)
 with threadpool_limits(limits=2):X=io.mmread(rd('matrix.mtx.gz')).tocsr()
 archive.close();assert X.shape==(len(ft),len(bar));umi=np.asarray(X.sum(0)).ravel();ng=np.asarray((X>0).sum(0)).ravel();mt=np.asarray(X[np.char.startswith(genes,'MT-')].sum(0)).ravel()/np.maximum(umi,1)
 if dataset=='PMGC2026':
  cells=pd.read_csv(O/'processed/PMGC2026_scRNA'/(sample+'_cell_source_QC.tsv.gz'),sep='\t');assert np.array_equal(cells.barcode.to_numpy(),bar);keep=cells.primary_QC.to_numpy(dtype=bool);labels=cells.group.to_numpy(dtype=str)
 else:
  counts={g:X[np.where(genes==g)[0][0]].toarray().ravel() for g in set(sum(panels.values(),[]))}
  def det(gs):return sum((counts[g]>0).astype(int) for g in gs)
  imm=det(panels['immune'])>=2;mes=(det(panels['mesothelial_specific'])>=1)&(det(panels['mesothelial_support'])>=1)&~imm;ep=(counts['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)&(det(panels['mesothelial_specific'])==0)&~imm;fib=(det(panels['fibroblast'])>=3)&~imm&~mes&~ep;end=(det(panels['endothelial'])>=2)&~imm&~mes&~ep;mur=(det(panels['mural'])>=2)&~imm&~mes&~ep;over=(fib.astype(int)+end.astype(int)+mur.astype(int))>1;groups={'immune_like':imm,'strict_mesothelial':mes,'conservative_epithelial':ep,'fibroblast_like':fib&~over,'endothelial_like':end&~over,'mural_like':mur&~over,'overlapping_stromal':over};groups['unclassified']=~np.logical_or.reduce(list(groups.values()));assert np.all(np.stack(list(groups.values())).sum(0)==1);labels=np.select(list(groups.values()),list(groups.keys()),default='ERROR');q=tissue['qc_primary'];keep=(umi>=q['UMI_min'])&(ng>=q['genes_min'])&(mt<=q['mitochondrial_fraction_max'])
 rng=np.random.default_rng(spec['sampling_seed']+case_index);selected=[];coverage={}
 for group,cap in spec['sampling_caps'].items():
  ids=np.flatnonzero(keep&(labels==group));take=ids if cap is None or len(ids)<=cap else rng.choice(ids,size=cap,replace=False);selected.extend(take.tolist());coverage[group]={'available':len(ids),'selected':len(take)}
 selected=np.array(sorted(selected));assert coverage['strict_mesothelial']['selected']>=20;assert all(coverage[g]['selected']>=20 for g in spec['primary_reference']);gene_keep=~np.isin(genes,spec['gene_exclusions']);counts=X[gene_keep][:,selected].T.tocsr().astype(np.float32);norm=sparse.diags(10000/umi[selected])@counts;norm=norm.tocsr();norm.data=np.log1p(norm.data).astype(np.float32);obs=pd.DataFrame({'group':labels[selected],'dataset':dataset,'sample':sample,'original_total_UMI':umi[selected],'original_detected_genes':ng[selected],'mitochondrial_fraction':mt[selected],'high_complexity_gt7000genes':ng[selected]>7000},index=pd.Index(bar[selected],name='barcode'));var=pd.DataFrame({'gene_symbol':genes[gene_keep]},index=pd.Index(ft[0].to_numpy(dtype=str)[gene_keep],name='gene_id'));a=ad.AnnData(norm,obs=obs,var=var);a.layers['counts']=counts;a.uns['normalization']='log1p(full-cell CP10k); target/gate genes excluded from CNVfeatures';p=out/(dataset+'_'+sample+'.h5ad');assert not p.exists();a.write_h5ad(p,compression='gzip');record={'dataset':dataset,'sample':sample,'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'cells':a.n_obs,'genes_before_coordinate_filter':a.n_vars,'sampling':coverage,'CGNL1_excluded':'CGNL1' not in a.var.gene_symbol.values,'all_gate_genes_excluded':not set(a.var.gene_symbol)&set(spec['gene_exclusions'])};assert record['CGNL1_excluded'] and record['all_gate_genes_excluded'];records.append(record);print('CNV_INPUT',json.dumps(record),flush=True)
(O/'metadata/CNV_INPUTS_READY_166.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cases':records,'spec_sha256':hashlib.sha256((O/'CNV_IDENTITY_SPECIFICATION_165.json').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
