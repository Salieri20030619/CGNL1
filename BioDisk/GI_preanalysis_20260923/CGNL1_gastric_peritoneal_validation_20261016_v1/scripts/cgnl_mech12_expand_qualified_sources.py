from pathlib import Path
import urllib.request,urllib.parse,requests,concurrent.futures,json,gzip,hashlib
O=Path.cwd()/'mechanism_investigation_v2'
queries={'mesothelial_gastric':'mesothelial AND (gastric OR stomach)','mesothelial_perturbation':'mesothelial AND (TGF OR cancer OR conditioned OR metastasis)','cgnl_exact':'CGNL1[Title] OR paracingulin[Title]','gastric_peritoneal_sc':'gastric AND peritoneal AND single cell'}
def search(k,q):
 r=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'gds','term':q,'retmode':'json','retmax':60},timeout=(15,60));p=O/'metadata'/('search_'+k+'.json');p.write_bytes(r.content);r.raise_for_status();d=r.json();print(k,d.get('esearchresult',{}).get('idlist'),flush=True);return {'key':k,'query':q,'status':r.status_code,'data':d}
out=[]
for k,q in queries.items():
 try:out.append(search(k,q))
 except Exception as e:out.append({'key':k,'error':str(e)})
u='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228955/suppl/GSE228955_Normalized_counts_cpm.txt.gz';r={'url':u}
try:
 with urllib.request.urlopen(u,timeout=45) as a:b=a.read()
 gzip.decompress(b);p=O/'raw/GSE228955_counts_ftp.gz';p.write_bytes(b);r.update(success=True,bytes=len(b),sha256=hashlib.sha256(b).hexdigest());print('FTP',r,flush=True)
except Exception as e:r.update(success=False,error=str(e));print('FTP',r,flush=True)
out.append(r);(O/'metadata/ACQUISITION_12.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
