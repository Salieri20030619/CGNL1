import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import pathlib,tarfile,gzip,io,json,re,hashlib,numpy as np,pandas as pd
from scipy import io as sio,sparse,stats
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2';O=R/'04_results'/'independent_scRNA';O.mkdir(parents=True,exist_ok=True)
lineage={'myeloid_core':['LST1','TYROBP','FCER1G','AIF1','CSF1R'],'macrophage_mature':['C1QA','C1QB','C1QC','CD68','MERTK','MS4A7','CD163','FCGR1A'],'T':['CD3D','CD3E','TRAC'],'B':['MS4A1','CD79A','CD79B'],'epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7'],'fibroblast':['COL1A1','COL1A2','DCN','LUM'],'endothelial':['PECAM1','VWF','KDR','CLDN5'],'mesothelial':['MSLN','CALB2','WT1','UPK3B']}
markers=['FABP4','CD36','LPL','FABP5','SLC27A1','SLC27A2','SLC27A4','LDLR','ADIPOQ','PLIN1','CIDEC','CIDEA','PPARG','TREM2','APOE','APOC1','LIPA','PLIN2','GPNMB','MARCO','MSR1','C1QA','SPP1','DAB2','VSIG4','FOLR2','MRC1','CD163','CD68','LST1','PTPRC','EPCAM','KRT19','COL1A1','DCN','PECAM1','MSLN','CALB2']
# Main validation endpoints fixed before external expression read; none used for gating.
endpoints=['CD36','TREM2','FABP4'];assert not set(endpoints)&set(sum(lineage.values(),[]))
(R/'00_docs'/'external_gating_plan.json').write_text(json.dumps({'lineage':lineage,'endpoints':endpoints,'scope':'Conservative marker-defined enriched compartments; not definitive single-cell atlas annotation. Gating excludes endpoint genes. Patient paired comparison is within same fluid sample; not metastasis vs healthy.'},indent=2))
soft=gzip.open(B/'01_metadata'/'GSE228598.soft.gz','rt').read();sampletypes={}
for z in soft.split('^SAMPLE = ')[1:]:
 acc=z.splitlines()[0];m=re.search(r'!Sample_characteristics_ch1 = sample type: (.+)',z);sampletypes[acc]=m.group(1).strip() if m else 'unknown'
archive=R/'02_raw'/'GSE228598_RAW.tar';out=[];qc=[];cellpaths=[]
with tarfile.open(archive) as tar:
 names=tar.getnames();prefixes=sorted(n[:-len('_matrix.mtx.gz')] for n in names if n.endswith('_matrix.mtx.gz'))
 for prefix in prefixes:
  acc=prefix.split('_')[0];patient=prefix.split('_')[1]
  def read(n):return gzip.GzipFile(fileobj=tar.extractfile(n))
  ft=pd.read_csv(read(prefix+'_features.tsv.gz'),sep='\t',header=None);genes=ft.iloc[:,1].astype(str).values;bar=pd.read_csv(read(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0].values
  with threadpool_limits(limits=2):X=sio.mmread(read(prefix+'_matrix.mtx.gz')).tocsr()
  assert X.shape==(len(genes),len(bar));assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data));umi=np.asarray(X.sum(axis=0)).ravel();ng=np.asarray((X>0).sum(axis=0)).ravel();mt=np.asarray(X[np.char.startswith(genes.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(umi,1)*100;keep=(umi>=500)&(ng>=200)&(ng<=7000)&(mt<=25)
  X=X[:,keep];bc=bar[keep];nn=X.shape[1]
  def detect(gs):return np.asarray((X[np.isin(genes,gs)]>0).sum(axis=0)).ravel()
  k={n:detect(gs) for n,gs in lineage.items()};my=k['myeloid_core']>=2;ly=(k['T']>=2)|(k['B']>=2);ep=k['epithelial']>=2
  label=np.repeat('Other_unresolved',nn).astype(object);label[ly&~my&~ep]='Lymphoid';label[my&~ly&~ep]='Myeloid_other';label[my&(k['macrophage_mature']>=2)&~ly&~ep]='Macrophage_enriched';label[ep&~my&~ly]='Epithelial_like';label[(k['fibroblast']>=2)&~ep&~my&~ly]='Fibroblast_like';label[(k['endothelial']>=2)&~ep&~my&~ly]='Endothelial_like';label[(k['mesothelial']>=2)&~my&~ly]='Mesothelial_like';label[(my&ly)|(my&ep)|(ly&ep)]='Mixed_lineage'
  cm=pd.DataFrame({'barcode':bc,'sample':prefix,'patient':patient,'sample_type':sampletypes[acc],'label':label,'umi':umi[keep],'n_genes':ng[keep],'mt_percent':mt[keep]});cm.to_csv(R/'03_processed'/(prefix+'_gating.tsv.gz'),sep='\t',index=False)
  for lab in np.unique(label):
   ix=np.flatnonzero(label==lab);total=float(X[:,ix].sum())
   for g in markers:
    gi=np.flatnonzero(genes==g)
    if len(gi):
     vals=np.asarray(X[gi][:,ix].sum(axis=0)).ravel();out.append({'patient':patient,'sample':prefix,'sample_type':sampletypes[acc],'celltype':lab,'gene':g,'n_cells':len(ix),'n_detected':int((vals>0).sum()),'fraction_detected':float((vals>0).mean()),'CPM':float(vals.sum()/max(total,1)*1e6)})
  qc.append({'patient':patient,'sample':prefix,'sample_type':sampletypes[acc],'input_cells':len(bar),'qc_cells':nn,'macrophage_enriched':int((label=='Macrophage_enriched').sum()),'myeloid_other':int((label=='Myeloid_other').sum()),'lymphoid':int((label=='Lymphoid').sum()),'epithelial_like':int((label=='Epithelial_like').sum()),'mixed_lineage':int((label=='Mixed_lineage').sum()),'genes':len(genes)})
  print(json.dumps(qc[-1]),flush=True)
a=pd.DataFrame(out);a['log2CPM']=np.log2(a.CPM+1);a.to_csv(O/'marker_expression_patient_compartment.tsv',sep='\t',index=False);q=pd.DataFrame(qc);q.to_csv(O/'sample_qc.tsv',sep='\t',index=False)
eligible=a[a.n_cells>=20];cmp=[]
for g in markers:
 z=eligible[(eligible.gene==g)&eligible.celltype.isin(['Macrophage_enriched','Lymphoid'])].pivot(index='patient',columns='celltype',values='log2CPM').dropna()
 if set(['Macrophage_enriched','Lymphoid']).issubset(z.columns):
  d=z.Macrophage_enriched-z.Lymphoid;n=len(d);pval=float(stats.wilcoxon(d,method='auto').pvalue) if n>=10 and (d!=0).any() and g not in set(sum(lineage.values(),[])) else np.nan
  cmp.append({'gene':g,'n_pairs':n,'median_log2CPM_difference':float(d.median()),'positive_pairs':int((d>0).sum()),'negative_pairs':int((d<0).sum()),'P':pval})
c=pd.DataFrame(cmp);c['FDR_BH_all_marker_family']=np.nan;ok=c.P.notna();c.loc[ok,'FDR_BH_all_marker_family']=stats.false_discovery_control(c.loc[ok,'P'].values,method='bh');c.to_csv(O/'within_patient_macrophage_vs_lymphoid.tsv',sep='\t',index=False)
end=c[c.gene.isin(endpoints)].copy();end['FDR_BH_prespecified_3_endpoints']=stats.false_discovery_control(end.P.values,method='bh');end.to_csv(O/'prespecified_endpoint_validation.tsv',sep='\t',index=False)
sumrows=[]
for g in endpoints:
 z=eligible[(eligible.celltype=='Macrophage_enriched')&(eligible.gene==g)];sumrows.append({'gene':g,'eligible_patients':int(z.patient.nunique()),'detected_in_patients':int(z[z.n_detected>0].patient.nunique()),'median_cell_detection_fraction':float(z.fraction_detected.median()),'median_log2CPM':float(z.log2CPM.median())})
s={'patients':int(q.patient.nunique()),'input_cells':int(q.input_cells.sum()),'qc_cells':int(q.qc_cells.sum()),'macrophage_enriched_cells':int(q.macrophage_enriched.sum()),'macrophage_eligible_patients':int((q.macrophage_enriched>=20).sum()),'endpoint_detection':sumrows,'scope':'Independent localization only. No primary tissue controls; ascites/washings not treated as metastatic vs nonmetastatic groups.'};(O/'summary.json').write_text(json.dumps(s,indent=2));(O/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'archive_readable':True}));print(json.dumps(s,indent=2));print(end.to_string(index=False))
