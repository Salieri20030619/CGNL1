import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import pathlib,json,gzip,io,re,hashlib,sys,itertools
import numpy as np,pandas as pd,h5py
from scipy import sparse,io as sio,stats
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');O=R/'04_results';F=O/'figures';F.mkdir(exist_ok=True)
def savefig(name):
 plt.tight_layout();plt.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');plt.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close()
def write(d,p):d.to_csv(p,sep='\t',index=False)
summary={}
# CRC effects and diagnostic panels: manually defined descriptive marker panels, no enrichment test.
C=O/'CRC';eff=pd.read_csv(C/'descriptive_gene_effects.tsv',sep='\t');q=pd.read_csv(C/'sample_qc.tsv',sep='\t')
panels={'Cell_cycle':['MKI67','PCNA','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','TOP2A','CCNB1','CDK1'], 'Translation':['RPL17','EEF1G','EEF2','EIF4A1','EIF4G1','RPS6','RPLP0','RPS3','RPL11'], 'Stress_response':['ATF3','ATF4','DDIT3','DDIT4','HSPA5','HMOX1','NQO1','TXNIP'], 'WNT_stem':['LGR5','ASCL2','AXIN2','OLFM4','MYC','SMOC2'], 'YAP_regeneration':['YAP1','WWTR1','CTGF','CYR61','ANKRD1','CLDN4','TACSTD2'], 'Lipid_synthesis':['FASN','ACACA','SCD','ACLY','SREBF1','ELOVL6','FADS2'], 'FAO':['CPT1A','CPT2','ACADM','ACADL','HADHA','HADHB','ACOX1'], 'Adhesion':['ITGA5','ITGB1','ITGAV','ITGB3','CDH1','CD44','PTK2','PXN','VCL']}
(R/'00_docs'/'descriptive_marker_panels.json').write_text(json.dumps({'scope':'Manually defined diagnostic panels; not validated signatures, not formal enrichment, no clinical score. Defined before reading gene-level effects.','panels':panels},indent=2))
rows=[]
for p,gs in panels.items():
 z=eff[eff.symbol.isin(gs)]
 for ct in ['DTP_vs_CTRL','Regrowth_vs_DTP','Regrowth_vs_CTRL','Regrowth_vs_DMSO_Barcode','Regrowth_vs_Saline']:
  rows.append({'panel':p,'contrast':ct,'genes_measured':len(z),'genes_planned':len(gs),'median_log2CPM_difference':z[ct].median(),'genes_positive':int((z[ct]>0).sum()),'genes_negative':int((z[ct]<0).sum()),'genes':';'.join(z.symbol)})
pa=pd.DataFrame(rows);write(pa,C/'diagnostic_panel_effects.tsv');write(eff[eff.symbol.isin(set(sum(panels.values(),[])))],C/'diagnostic_panel_genes.tsv')
pca=pd.read_csv(C/'pca.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,5))
for g,z in pca.groupby('group'):ax.scatter(z.PC1,z.PC2,label=g,s=40)
ax.set(xlabel='PC1',ylabel='PC2',title='CRC: 25 PDX libraries; one donor for main contrasts');ax.legend(fontsize=7);savefig('CRC_01_PCA')
piv=pa.pivot(index='panel',columns='contrast',values='median_log2CPM_difference')[['DTP_vs_CTRL','Regrowth_vs_DTP','Regrowth_vs_CTRL']]
fig,ax=plt.subplots(figsize=(7,5));im=ax.imshow(piv.values,cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(len(piv)),piv.index);ax.set_xticks(range(3),['DTP - control','Regrowth - DTP','Regrowth - control'],rotation=20,ha='right');fig.colorbar(im,ax=ax,label='Median log2 CPM difference');ax.set_title('CRC descriptive panels: no patient-level inference');savefig('CRC_02_panel_effects')
summary['CRC']=dict(pd.read_csv(C/'summary.tsv',sep='\t').values)
# Multiome structural and count QC. No fragments/TSS/FRiP; not a regulatory validation.
M=O/'CRC_multiome';M.mkdir(exist_ok=True);multi=[]
for p in sorted((R/'02_raw').glob('GSE294559_*.h5')):
 with h5py.File(p) as h:
  m=h['matrix'];x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));types=m['features']['feature_type'][:].astype(str);names=m['features']['name'][:].astype(str)
  xr=x[types=='Gene Expression'];xp=x[types=='Peaks'];gene_names=names[types=='Gene Expression'];umi=np.asarray(xr.sum(axis=0)).ravel();ng=np.asarray((xr>0).sum(axis=0)).ravel();at=np.asarray(xp.sum(axis=0)).ravel();npk=np.asarray((xp>0).sum(axis=0)).ravel();mt=np.asarray(xr[np.char.startswith(gene_names,'MT-')].sum(axis=0)).ravel()/np.maximum(umi,1)*100
  met=pd.DataFrame({'barcode':m['barcodes'][:].astype(str),'rna_umi':umi,'rna_genes':ng,'atac_peak_counts':at,'atac_detected_peaks':npk,'mt_percent':mt});met['basic_joint_pass']=(ng>=200)&(umi>=500)&(mt<=25)&(npk>=500)
  write(met,M/(p.stem+'_barcode_qc.tsv'))
  multi.append({'file':p.name,'patient':'P001','type':'tumor_organoid' if '_TO' in p.name else 'tumor_tissue','n_barcodes':x.shape[1],'gene_features':xr.shape[0],'peak_features':xp.shape[0],'basic_joint_pass':int(met.basic_joint_pass.sum()),'median_rna_genes':float(np.median(ng)),'median_atac_peaks':float(np.median(npk)),'treatment_timecourse':False})
write(pd.DataFrame(multi),M/'multiome_structure_qc.tsv');summary['CRC_multiome']=multi
# GC paired scRNA: verify supplied annotation, sparsity and marker expression without inventing malignant cells.
G=O/'GC';G.mkdir(exist_ok=True);meta=pd.read_csv(R/'02_raw'/'GSE234129_meta.tsv.gz',sep='\t');bar=pd.read_csv(R/'02_raw'/'GSE234129_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
features=pd.read_csv(R/'02_raw'/'GSE234129_features.tsv.gz',sep='\t',header=None);symbols=features.iloc[:,1 if features.shape[1]>1 else 0].astype(str).values
with threadpool_limits(limits=2):x=sio.mmread(R/'02_raw'/'GSE234129_count_matrix.mtx.gz').tocsr()
assert x.shape==(len(symbols),len(bar)),(x.shape,len(symbols),len(bar));assert not bar.duplicated().any();assert set(bar)==set(meta.cell_barcodes);meta=meta.set_index('cell_barcodes').loc[bar].reset_index()
assert np.all(x.data>=0) and np.all(x.data==np.round(x.data))
umi=np.asarray(x.sum(axis=0)).ravel();ng=np.asarray((x>0).sum(axis=0)).ravel();mt=np.asarray(x[np.char.startswith(symbols.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(umi,1)*100
meta['n_umi']=umi;meta['n_genes']=ng;meta['mt_percent']=mt;meta['basic_pass']=(ng>=200)&(umi>=500)&(mt<=25)
write(meta,G/'GSE234129_cell_qc.tsv');ct=pd.crosstab(meta['sample'],meta['celltype']);ct.to_csv(G/'GSE234129_sample_celltype_counts.tsv',sep='\t')
markers=['EPCAM','KRT8','KRT18','KRT19','PTPRC','CD3D','LST1','COL1A1','MSLN','CALB2','WT1'];mr=[]
for s,ix in meta.groupby('sample').groups.items():
 for g in markers:
  gi=np.flatnonzero(symbols==g)
  if len(gi):a=np.asarray(x[gi][:,list(ix)].sum(axis=0)).ravel();mr.append({'sample':s,'gene':g,'n_cells':len(ix),'fraction_detected':float(np.mean(a>0)),'mean_count':float(a.mean())})
write(pd.DataFrame(mr),G/'GSE234129_marker_detection.tsv')
annot_ep=meta.celltype.str.contains('epithel|malignan|tumou?r|cancer',case=False,regex=True)
summary['GC_scRNA']={'genes':x.shape[0],'cells':x.shape[1],'patients':int(meta.patient.nunique()),'samples':int(meta['sample'].nunique()),'annotated_epithelial_or_malignant_cells':int(annot_ep.sum()),'basic_qc_pass':int(meta.basic_pass.sum()),'primary_ascites_pairs_in_design':4,'usable_malignant_pairs_from_supplied_annotation':0,'mesothelial_cells':int(meta.celltype.str.startswith('Mesothelial').sum())}
fig,ax=plt.subplots(figsize=(8,5));gct=meta.groupby('sample').size();ax.barh(gct.index,gct.values,color='#4C78A8');ax.set(xlabel='Cells in public matrix',title='GC paired dataset: immune/stromal annotations; no malignant labels');savefig('GC_01_cells_by_sample')
# GC 3-pair microarray. Gene symbols absent in deposited platform; keep probe-level results.
text=gzip.open(R/'02_raw'/'GSE200697_family.soft.gz','rt').read().splitlines();samples={};annotation=None;cur=None;i=0
while i<len(text):
 l=text[i]
 if l.startswith('^SAMPLE = '):cur=l.split(' = ')[1];samples[cur]={'meta':{},'table':None}
 elif cur is not None and l.startswith('!Sample_') and ' = ' in l:
  k,v=l.split(' = ',1);samples[cur]['meta'].setdefault(k,[]).append(v)
 elif l=='!platform_table_begin' or l=='!sample_table_begin':
  end='!platform_table_end' if l=='!platform_table_begin' else '!sample_table_end';j=i+1
  while text[j]!=end:j+=1
  df=pd.read_csv(io.StringIO('\n'.join(text[i+1:j])),sep='\t',low_memory=False)
  if end=='!platform_table_end':annotation=df
  else:samples[cur]['table']=df
  i=j
 i+=1
ann=annotation.set_index('ID');cols=[];design=[]
for acc,s in samples.items():
 t=s['table'];cols.append(t.groupby('ID_REF')['VALUE'].median().rename(acc));title=s['meta']['!Sample_title'][0];pair=re.search(r'Tissue(\d+)',title).group(1);condition='PM' if 'peritoneal' in title else 'Primary';design.append({'sample':acc,'pair':pair,'condition':condition,'title':title,'clinical':';'.join(s['meta'].get('!Sample_characteristics_ch1',[]))})
a=pd.concat(cols,axis=1);des=pd.DataFrame(design);write(des,G/'GSE200697_pair_design.tsv');a.to_csv(R/'03_processed'/'GC_microarray_processed_probe_values.tsv.gz',sep='\t',compression='gzip');write(annotation,G/'GSE200697_platform_annotation.tsv')
# Values are distributed on log-like scale; retain supplied scale and do not mislabel as fold changes.
diffs=[]
for pair,z in des.groupby('pair'):
 pm=z[z.condition=='PM']['sample'].iloc[0];pr=z[z.condition=='Primary']['sample'].iloc[0];diffs.append((a[pm]-a[pr]).rename('pair_'+pair))
delta=pd.concat(diffs,axis=1);res=delta.copy();res['mean_difference_supplied_scale']=delta.mean(axis=1);res['median_difference_supplied_scale']=delta.median(axis=1);res['positive_pairs']=(delta>0).sum(axis=1);res['negative_pairs']=(delta<0).sum(axis=1);res['transcript_type']=ann['TRANSCRIPT_TYPE'].reindex(res.index);res['annotation_available']=res.index.isin(ann.index)
valid=res.annotation_available&res.transcript_type.notna();res=res[valid].copy();dv=delta.loc[res.index].values;obs=np.abs(dv.mean(axis=1));ex=np.zeros(len(res))
for signs in itertools.product([-1,1],repeat=3):ex+=(np.abs((dv*np.array(signs)).mean(axis=1))>=obs-1e-12)
res['P_exact_signflip']=ex/8;res['FDR_BH']=np.nan
for typ,ix in res.groupby('transcript_type').groups.items():res.loc[ix,'FDR_BH']=stats.false_discovery_control(res.loc[ix,'P_exact_signflip'].values,method='bh')
res.index.name='probe';res.to_csv(G/'GSE200697_probe_effects.tsv',sep='\t');a.describe().to_csv(G/'GSE200697_sample_distributions.tsv',sep='\t')
summary['GC_microarray']={'unique_probes_in_matrix':len(a),'annotated_probes_analyzed':len(res),'samples':len(des),'independent_pairs':3,'minimum_exact_two_sided_P':.25,'FDR_lt_005':int((res.FDR_BH<.05).sum()),'direction_consistent_abs_mean_ge1':int((((res.positive_pairs==3)|(res.negative_pairs==3))&(res.mean_difference_supplied_scale.abs()>=1)).sum()),'platform_columns':annotation.columns.tolist(),'gene_symbols_available':False,'expression_scale':'submitted normalized values, not independently confirmed log base'}
# Protein resource QC. Five letter-coded groups lack meaning and biological replicate map.
p=pd.read_csv(R/'02_raw'/'PXD064711_protein.raw.intensity.txt',sep='\t',encoding='utf-16',low_memory=False);intcols=[c for c in p if re.fullmatch(r'[ABCDE][1-5]|ref[12]',c)];v=p[intcols].apply(pd.to_numeric,errors='coerce');qc=pd.DataFrame({'sample':intcols,'positive_proteins':(v>0).sum().values,'missing_or_zero':(v.isna()|(v<=0)).sum().values,'median_positive_intensity':v.where(v>0).median().values});write(qc,G/'PXD064711_sample_qc.tsv');write(p[['ID','ProteinGroup','GeneName']+intcols],R/'03_processed'/'GC_protein_quantification.tsv.gz')
summary['GC_protein']={'protein_rows':len(p),'quantitative_channels':len(intcols),'biological_label_channels':len([c for c in intcols if not c.startswith('ref')]),'reference_channels':2,'groups':['A','B','C','D','E'],'group_meanings_verified':False,'biological_replicate_identity_verified':False,'gene_names_nonmissing':int(p.GeneName.notna().sum()),'differential_testing':'not performed without group mapping'}
fig,ax=plt.subplots(figsize=(8,4));ax.bar(qc['sample'],qc.positive_proteins,color='#59A14F');ax.set(xlabel='Submitted channel labels; biological groups unknown',ylabel='Proteins with positive intensity',title='GC protein table: technical usability only');savefig('GC_02_protein_coverage')
(O/'preanalysis_summary.json').write_text(json.dumps(summary,indent=2));(R/'logs'/'python_environment.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'h5py':h5py.__version__},indent=2));(O/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'analyses':['CRC descriptive','CRC multiome structural QC','GC scRNA QC','GC paired probe effects','GC protein QC']}));print(json.dumps(summary,indent=2))
