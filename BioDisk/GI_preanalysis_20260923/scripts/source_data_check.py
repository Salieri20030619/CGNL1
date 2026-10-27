import pathlib,re,json,pandas as pd
from bs4 import BeautifulSoup
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
p=r/'01_metadata'/'GC2026.html';s=BeautifulSoup(p.read_text(),'html.parser')
for sec in s.find_all(['section','div']):
 if sec.get('id','') in ['data-availability','data-availability-section']:print('GC2026_DATA',sec.get_text(' ',strip=True))
t=s.get_text(' ',strip=True)
for key in ['Data availability','GSE','EGAS']:
 idx=t.find(key);print(key,t[idx:idx+1500])
print('GC2026_SUPPL',[(a.get_text(' ',strip=True),a.get('href')) for a in s.find_all('a') if 'MOESM' in a.get('href','')])
p=r/'01_metadata'/'GC_ascites_source.xlsx';x=pd.ExcelFile(p);print('ASCITES_SOURCE_SHEETS',x.sheet_names)
for sh in x.sheet_names:
 d=pd.read_excel(p,sheet_name=sh);print(sh,d.shape,d.head(2).to_string(index=False)[:1000])
p=r/'01_metadata'/'iprox_project.html';t=p.read_text();print('IPROX_INLINE',t[-12000:][:6000])
