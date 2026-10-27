import pathlib,subprocess,concurrent.futures,xml.etree.ElementTree as E,json,gzip,csv,collections
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');M=O/'metadata'
r=E.parse(M/'GSE225857_article.xml').getroot()
for p in r.iter('p'):
 t=''.join(p.itertext())
 if any(x in t for x in ['six patients','6 patients','Supplementary Table 1','malignant','copy number','C1','s1231']):print(t[:1800])
for x in r.iter('supplementary-material'):print(E.tostring(x,encoding='unicode')[:1200])
with gzip.open(M/'nonimmune_meta.txt.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
print('CELLS',len(rows));print('PATIENT_ORGAN',collections.Counter((r['patients'],r['organs']) for r in rows));print('CLUSTERS',collections.Counter(r['cluster'] for r in rows))
jobs=[('visium_meta_after_qc.rds','https://raw.githubusercontent.com/yliuup/CRC_micromets_ST/main/Meta_data/visium_meta_after_qc.rds'),('visium_Readme.md','https://raw.githubusercontent.com/yliuup/CRC_micromets_ST/main/README.md'),('Figure1D_E.R','https://raw.githubusercontent.com/yliuup/CRC_micromets_ST/main/Figure1/Figure1D_E.R'),('spatial_table1.xlsx','https://ars.els-cdn.com/content/image/1-s2.0-S1535610826002965-mmc1.xlsx'),('spatial_table2.xlsx','https://ars.els-cdn.com/content/image/1-s2.0-S1535610826002965-mmc2.xlsx')]
for x in r.iter('supplementary-material'):
 for q in x.iter():
  u=q.attrib.get('{http://www.w3.org/1999/xlink}href','')
  if u.endswith(('.xlsx','.pdf','.docx')):jobs.append((u,'https://pmc.ncbi.nlm.nih.gov/articles/instance/10275599/bin/'+u))
def run(v):
 n,u=v;r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','90','-o',str(M/n),u],capture_output=True);return n,r.returncode,(M/n).stat().st_size if (M/n).exists() else 0
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:print(list(ex.map(run,jobs)))
