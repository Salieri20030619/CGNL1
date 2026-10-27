import sys,concurrent.futures
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/scripts')
from acquire import fetch
items=[('CRC_diapause.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8437243/fullTextXML'),('GC_2021.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8074162/fullTextXML'),('GC_228598.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11012629/fullTextXML'),('PXD064711_central.xml','https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD064711&outputMode=XML'),('GC_ascites_source.xlsx','https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-023-36310-9/MediaObjects/41467_2023_36310_MOESM4_ESM.xlsx'),('GSE145356_counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145356/suppl/GSE145356_merged_HTSeq_readCounts_V2.txt.gz')]
for a in ['GSE228598','GSE163558','GSE234129','GSE200697']:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}'
 items += [(a+'.soft.gz',base+'/soft/'+a+'_family.soft.gz'),(a+'.files.html',base+'/suppl/')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(fetch,items))
