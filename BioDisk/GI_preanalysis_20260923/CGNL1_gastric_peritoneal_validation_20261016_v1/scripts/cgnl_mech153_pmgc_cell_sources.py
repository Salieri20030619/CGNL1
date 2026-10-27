from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,zipfile,gzip,hashlib,datetime,pandas as pd,numpy as np
from scipy import io
from threadpoolctl import threadpool_limits
O=Path.cwd()/'mechanism_investigation_v2';q=json.loads((O/'metadata/PMGC2026_IMPORT_Human_scRNA.zip.json').read_text());assert q['success'];p=Path(q['file']);assert hashlib.md5(p.read_bytes()).hexdigest()==q['md5'];spec=json.loads((O/'PMGC2026_QC_SPECIFICATION_145.json').read_text());panels=json.loads((O/'TISSUE_ANALYSIS_SPECIFICATION.json').read_text())['source_panels'];needed=set(sum(panels.values(),[])+['CGNL1']);dest=O/'processed/PMGC2026_scRNA';dest.mkdir(exist_ok=True);rows=[];qc=[];features=[]
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None;sample_names=sorted({n.split('/')[1] for n in z.namelist() if n.endswith('/matrix.mtx.gz')});assert sample_names==['PM1','PM2','PM3','PM4']
 for sample in sample_names:
  prefix='Human_scRNA/'+sample+'/'
  def rd(name):return gzip.GzipFile(fileobj=z.open(prefix+name))
  ft=pd.read_csv(rd('features.tsv.gz'),sep='\t',header=None);bar=pd.read_csv(rd('barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str);genes=ft[1].astype(str).to_numpy();assert bar.is_unique and ft[0].is_unique;assert ft.loc[ft[1]=='CGNL1',0].tolist()==['ENSG00000128849'];assert ft[0].str.startswith('ENSG').all()
  with threadpool_limits(limits=2):X=io.mmread(rd('matrix.mtx.gz')).tocsr()
  assert X.shape==(len(genes),len(bar)) and (X.data>=0).all() and np.equal(X.data,np.floor(X.data)).all();umi=np.asarray(X.sum(0)).ravel();ng=np.asarray((X>0).sum(0)).ravel();mito=np.asarray(X[np.char.startswith(genes,'MT-')].sum(0)).ravel()/np.maximum(umi,1);hb=spec['primary_scRNA_QC']['erythrocyte_gene_set'];ery=np.asarray(X[np.isin(genes,hb)].sum(0)).ravel()/np.maximum(umi,1);counts={}
  for g in needed:
   ix=np.where(genes==g)[0];assert len(ix)==1,(sample,g,len(ix));counts[g]=X[ix[0]].toarray().ravel()
  features.append({'sample':sample,'genes':len(genes),'cells':len(bar),'target_present':True,'missing_erythrocyte_genes':[g for g in hb if g not in genes],'UMI_sum':int(umi.sum()),'integer_counts':True})
  def det(gs):return sum((counts[g]>0).astype(int) for g in gs)
  immune=det(panels['immune'])>=2;meso=(det(panels['mesothelial_specific'])>=1)&(det(panels['mesothelial_support'])>=1)&~immune;epi=(counts['EPCAM']>0)&(det(['KRT8','KRT18','KRT19'])>=2)&(det(panels['mesothelial_specific'])==0)&~immune;fib=(det(panels['fibroblast'])>=3)&~immune&~meso&~epi;end=(det(panels['endothelial'])>=2)&~immune&~meso&~epi;mural=(det(panels['mural'])>=2)&~immune&~meso&~epi;overlap=(fib.astype(int)+end.astype(int)+mural.astype(int))>1;groups={'immune_like':immune,'strict_mesothelial':meso,'conservative_epithelial':epi,'fibroblast_like':fib&~overlap,'endothelial_like':end&~overlap,'mural_like':mural&~overlap,'overlapping_stromal':overlap};groups['unclassified']=~np.logical_or.reduce(list(groups.values()));assert np.all(np.stack(list(groups.values())).sum(0)==1)
  c=spec['primary_scRNA_QC'];primary=(ng>=c['genes_min_inclusive'])&(ng<=c['genes_max_inclusive'])&(umi<c['UMI_max_exclusive'])&(mito<c['mitochondrial_fraction_max_exclusive'])&(ery<c['erythrocyte_fraction_max_exclusive']);sensitivity=(umi>=500)&(ng>=200)&(ng<=7000)&(mito<=.25);cell=pd.DataFrame({'barcode':bar,'sample':sample,'total_UMI':umi,'detected_genes':ng,'mitochondrial_fraction':mito,'erythrocyte_fraction':ery,'primary_QC':primary,'sensitivity_QC':sensitivity,'group':np.select(list(groups.values()),list(groups.keys()),default='ERROR'),**counts});cell.to_csv(dest/(sample+'_cell_source_QC.tsv.gz'),sep='\t',index=False,encoding='utf-8')
  for scheme,keep in [('primary',primary),('sensitivity',sensitivity)]:
   qc.append({'sample':sample,'scheme':scheme,'raw_cells':len(bar),'retained_cells':int(keep.sum()),'CGNL1_UMI':int(counts['CGNL1'][keep].sum()),'doublet_filter_reproduced':False})
   for group,gmask in groups.items():
    mask=gmask&keep;n=int(mask.sum());tot=int(umi[mask].sum());tar=counts['CGNL1'][mask];rows.append({'sample':sample,'scheme':scheme,'group':group,'n_cells':n,'total_UMI':tot,'CGNL1_UMI':int(tar.sum()),'CGNL1_positive_cells':int((tar>0).sum()),'CGNL1_detection':float((tar>0).mean()) if n else np.nan,'CGNL1_ln1pCP10k':float(np.log1p(tar.sum()/tot*1e4)) if tot else np.nan,'coverage20':n>=20,'coverage50':n>=50})
  print('PMGC_SOURCE',sample,qc[-2:],flush=True)
pd.DataFrame(rows).to_csv(O/'results/PMGC2026_scRNA_source_profiles.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(qc).to_csv(O/'results/PMGC2026_scRNA_QC.tsv',sep='\t',index=False,encoding='utf-8');out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':features,'scRNA_donors_per_deposit':4,'no_inferential_P':True,'source_identity':'targetexcludedmarkerqualified; noDNA/protein confirmation','author_labels_provided':False,'independent_numeric_QA_pending':True};(O/'results/PMGC2026_scRNA_SOURCE_COMPLETE.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
