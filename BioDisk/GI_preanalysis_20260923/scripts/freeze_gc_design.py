import pathlib,h5py,pandas as pd,json
from bs4 import BeautifulSoup
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
with h5py.File(r/'02_raw'/'GSE294559_TP001.h5') as h:
 f=h['matrix']['features'];k=f['feature_type'][:].astype(str)=='Gene Expression';d=pd.DataFrame({'gene_id':f['id'][:].astype(str)[k],'symbol':f['name'][:].astype(str)[k]});d.to_csv(r/'01_metadata'/'human_gene_symbols_from_multiome.tsv',sep='\t',index=False)
m=pd.read_csv(r/'01_metadata'/'GSE314812_sample_design.tsv',sep='\t');m=m[m.site.isin(['GCPT','GCPM'])];m=m[m.donor_label!='GC05b'];ct=pd.crosstab(m.patient,m.site);ids=ct.index[(ct.GCPT==1)&(ct.GCPM==1)];m=m[m.patient.isin(ids)].copy();m['role']='primary_paired';m.to_csv(r/'01_metadata'/'GSE314812_paired_design.tsv',sep='\t',index=False);print('PAIRED',len(ids),len(m),ids.tolist());print('PRESERVATION_CROSS',pd.crosstab(m.site,m['preservation method']).to_string())
# Mutation table context and tissue purity considerations
s=BeautifulSoup((r/'01_metadata'/'GC2026.html').read_text(),'html.parser').get_text(' ',strip=True)
for key in ['Supplementary Data 3','Supplementary Data 4','Supplementary Data 5','Supplementary Data 6','purity']:
 idx=s.find(key);print(key,s[max(0,idx-400):idx+150])
print('CLINICAL_HEADERS',pd.read_excel(r/'01_metadata'/'GC2026_clinical.xlsx',header=2).columns.tolist())
