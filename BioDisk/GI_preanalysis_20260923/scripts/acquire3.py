import sys,concurrent.futures
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/scripts')
from acquire import fetch
items=[]
for a in ['GSE163171','GSE187650','GSE163035','GSE187512','GSE294559']:
 b=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}'
 items.extend([(a+'.soft.gz',b+'/soft/'+a+'_family.soft.gz'),(a+'.files.html',b+'/suppl/')])
items.extend([('GC2026.html','https://www.nature.com/articles/s42003-026-09906-2'),('CRC_regrowth_2025.html','https://www.nature.com/articles/s41420-025-02567-w')])
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(fetch,items))
