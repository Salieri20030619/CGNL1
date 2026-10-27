import pathlib,tarfile,gzip,re,json,pandas as pd
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');s=gzip.open(O/'metadata/GSE50760.soft.gz','rt').read();rows=[]
for c in s.split('^SAMPLE = ')[1:]:
 acc=c.splitlines()[0];title=re.search(r'!Sample_title = (.+)',c).group(1).strip();pt=re.search(r'AMC_(\d+)',title).group(1);tissue=re.search(r'!Sample_characteristics_ch1 = tissue: (.+)',c).group(1).strip();site='PT' if tissue=='primary colorectal cancer' else ('LM' if tissue=='metastatic colorectal cancer to the liver' else 'NC');assert site!='NC' or tissue=='normal-looking surrounding colonic epithelium';rows.append({'sample':acc,'patient':pt,'site':site,'title':title,'colon_rectum':'not provided in GEO'})
m=pd.DataFrame(rows);assert len(m)==54 and m.patient.nunique()==18;assert (m.groupby(['patient','site']).size()==1).all();m.to_csv(O/'metadata/bulk_design.tsv',sep='\t',index=False)
vals=[]
with tarfile.open(O/'raw/GSE50760_RAW.tar') as tar:
 for x in tar.getmembers():
  if not x.name.endswith('.txt.gz'):continue
  f=pd.read_csv(gzip.GzipFile(fileobj=tar.extractfile(x)),sep='\t');bad=f.iloc[:,0].astype(str).str.match(r'^\d{1,2}[-/](Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)',case=False);f=f.loc[~bad];assert f.iloc[:,0].is_unique;assert (f.iloc[:,1]>=0).all();v=f.set_index(f.columns[0]).iloc[:,0];v.name=x.name.split('_')[0];vals.append(v)
a=pd.concat(vals,axis=1);assert not a.isna().any().any();a.to_csv(O/'raw/bulk_FPKM.tsv.gz',sep='\t',index_label='gene');print('FPKM',a.shape);print(m.groupby('site').size().to_string());plan=json.loads((O/'PLAN.json').read_text());gs=plan['frozen_genes'];print('CANDIDATE_COVERAGE',sum(g.replace('.','-') in a.index for g in gs),'MISSING',[g for g in gs if g.replace('.','-') not in a.index]);print(m.head().to_string(index=False))
