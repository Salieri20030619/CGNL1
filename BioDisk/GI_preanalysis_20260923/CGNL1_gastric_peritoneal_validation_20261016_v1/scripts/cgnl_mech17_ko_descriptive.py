from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,hashlib,numpy as np,pandas as pd
O=Path.cwd()/'mechanism_investigation_v2';p=O/'raw/GSE228955_counts_ftp.gz'
panels={'target_paralog':['CGNL1','CGN'],'junction_components':['CLDN1','CLDN2','CLDN3','CLDN4','CLDN7','OCLN','TJP1','TJP2','TJP3','CDH1','CTNNB1','JUP','PLEKHA7'],'microtubule_anchor':['CAMSAP3','CAMSAP2','MAPRE1','MAPRE3','KIFC3'],'myosin_contractility':['MYH9','MYH10','MYL9','MYL12A','MYL12B','ROCK1','ROCK2'],'Rho_Rac_regulators':['ARHGEF2','TIAM1','RHOA','RAC1','CDC42','RACGAP1','SH3BP1'],'mesenchymal':['VIM','FN1','COL1A1','COL1A2','ACTA2','TAGLN'],'cycling':['MKI67','TOP2A','PCNA','UBE2C','CENPF']}
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scale':'author normalized CPM; analysis log2(CPM+1), not raw counts','contrasts':['CGNL1_KO-WT','CGN_KO-WT','double_KO-WT','CGNL1_KO-CGN_KO'],'panels':panels,'analysis':'all measured genes preserve group values and mean log2 difference; predefined machinery subset for interpretation','no_inferential_P':'one clone per genotype; independent culture repeat provenance not sufficiently resolved; no P/FDR computed this stage','limitations':['Canine kidney epithelial model, not human mesothelium or gastric carcinoma','RNA abundance cannot test protein localization, GTPase activity, contractile force or actual barrier','No genotype-specific rescue transcriptomes','No selecting and calling large-fold genes significant without valid P']}
sp=O/'KO_DESCRIPTIVE_SPECIFICATION.json';assert not sp.exists();sp.write_text(json.dumps(spec,indent=2),encoding='utf-8')
x=pd.read_csv(p,sep='\t',index_col=0);assert x.shape==(13672,12);assert x.index.is_unique;assert np.isfinite(x.values).all() and (x.values>=0).all()
groups={'CGNL1_KO':[c for c in x if c.startswith('CGNL1_KO')],'CGN_KO':[c for c in x if c.startswith('CGN_KO')],'double_KO':[c for c in x if c.startswith('CGN_CGNL1')],'WT':[c for c in x if c.startswith('WT')]};assert all(len(v)==3 for v in groups.values());assert sum(groups.values(),[]).__len__()==12
y=np.log2(x+1);a=x.copy();a.index.name='gene'
for g,cs in groups.items():a[g+'_mean_CPM']=x[cs].mean(axis=1);a[g+'_mean_log2CPM1']=y[cs].mean(axis=1);a[g+'_SD_log2CPM1']=y[cs].std(axis=1)
for g in ['CGNL1_KO','CGN_KO','double_KO']:
 a[g+'_minus_WT']=a[g+'_mean_log2CPM1']-a.WT_mean_log2CPM1
 a[g+'_range_above_WT']=y[groups[g]].min(axis=1)>y[groups['WT']].max(axis=1)
 a[g+'_range_below_WT']=y[groups[g]].max(axis=1)<y[groups['WT']].min(axis=1)
a['CGNL1_minus_CGN']=a.CGNL1_KO_mean_log2CPM1-a.CGN_KO_mean_log2CPM1
a.to_csv(O/'results/GSE228955_all_gene_descriptive.tsv.gz',sep='\t',encoding='utf-8')
rows=[]
for panel,genes in panels.items():
 for g in genes:
  if g in a.index:rows.append(dict(panel=panel,gene=g,**a.loc[g].to_dict()))
  else:rows.append({'panel':panel,'gene':g,'status':'not_in_processed_matrix'})
t=pd.DataFrame(rows);t.to_csv(O/'results/GSE228955_prespecified_machinery.tsv',sep='\t',index=False,encoding='utf-8')
cols=['panel','gene','WT_mean_CPM','CGNL1_KO_mean_CPM','CGNL1_KO_minus_WT','CGN_KO_minus_WT','double_KO_minus_WT']
print(t[cols].to_string(index=False));receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'groups':groups,'genes':len(a),'spec_sha256':hashlib.sha256(sp.read_bytes()).hexdigest(),'matrix_verified':True,'P_values_computed':False,'mechanism_established':False};(O/'results/KO_DESCRIPTIVE_COMPLETE.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
