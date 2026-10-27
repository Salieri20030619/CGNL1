import pathlib,gzip,re,pandas as pd
from bs4 import BeautifulSoup
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');records=[];cur=None
for l in gzip.open(r/'01_metadata'/'GSE314812.soft.gz','rt'):
 l=l.strip()
 if l.startswith('^SAMPLE = '):
  if cur:records.append(cur)
  cur={'accession':l.split(' = ')[1]}
 elif cur is not None:
  if l.startswith('!Sample_title = '):cur['title']=l.split(' = ')[1];cur['donor_label']=re.search(r'patient (\S+)',cur['title']).group(1);cur['patient']=re.sub('[ab]$','',cur['donor_label'])
  if l.startswith('!Sample_description = Library name: '):cur['sample']=l.split('Library name: ')[1]
  if l.startswith('!Sample_characteristics_ch1 = '):
   z=l.split(' = ')[1];k,v=z.split(': ',1);cur[k]=v
if cur:records.append(cur)
m=pd.DataFrame(records);m['site']=m['sample'].str.split('_').str[-1];m.to_csv(r/'01_metadata'/'GSE314812_sample_design.tsv',sep='\t',index=False)
x=pd.read_csv(r/'02_raw'/'GSE314812_counts.tsv.gz',sep='\t');print('COUNT_SHAPE',x.shape,'COLUMNS',x.columns.tolist());print(x.iloc[:3,:5].to_string(index=False));print('DESIGN',pd.crosstab(m.patient,m.site).to_string());print('PRESERVATION',pd.crosstab(m.site,m['preservation method']).to_string());print('EXTRA_COLS',set(x.columns)-set(m['sample']))
t=BeautifulSoup((r/'01_metadata'/'GC2026.html').read_text(),'html.parser').get_text(' ',strip=True)
for key in ['chemotherapy','GC05','FFPE','treatment-na','treatment']:
 positions=[z.start() for z in re.finditer(re.escape(key),t,flags=re.I)]
 print('ARTICLE',key,[t[max(0,i-140):i+500] for i in positions[:5]])
