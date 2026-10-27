from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2'
urls={
 'GSE308231_family.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308231&targ=all&form=text&view=full',
 'GSE308231_filelist.txt':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308231/suppl/filelist.txt',
 'HRA002712.html':'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA002712',
 'ascites_suppdata.xlsx':'https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-023-36310-9/MediaObjects/41467_2023_36310_MOESM4_ESM.xlsx',
 'ascites_sourcedata.xlsx':'https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-023-36310-9/MediaObjects/41467_2023_36310_MOESM5_ESM.xlsx',
 'GSE308231_paper.json':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:41360923&resultType=core&format=json',
 'cgnl_geo_search.json':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=CGNL1%20OR%20paracingulin&retmode=json&retmax=50',
}
def get(k,u):
 p=O/'metadata'/k;r={'name':k,'url':u,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  a=requests.get(u,timeout=(15,60));p.write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest());print(k,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=str(e);print(k,str(e),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:receipt=list(ex.map(lambda kv:get(*kv),urls.items()))
(O/'metadata/ACQUISITION_03.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
