from pathlib import Path
import json,datetime,hashlib,concurrent.futures,requests
O=Path.cwd();N=O/'mechanism_investigation_v2'
for d in ['metadata','literature','raw','results','scripts','logs','report']: (N/d).mkdir(parents=True,exist_ok=True)
def put(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
put(N/'BASELINE.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'authorization':'Continue sustained investigation toward evidence-supported mechanism mainline; download actual data; retain all records','prior_decision':json.loads((O/'DECISION.json').read_text(encoding='utf-8')),'prior_report_sha256':hashlib.sha256((O/'report/CGNL1_gastric_peritoneal_validation.html').read_bytes()).hexdigest(),'criteria':['Verified patient-level disease comparison within a qualified cell source','Independent evidence separating source abundance from within-source state','Functional evidence specific to CGNL1, with direction and rescue where available','Novelty evaluated against existing primary literature','Computational perturbation is model evidence, not experimental causation'],'competing_explanations':['normal mesothelial identity / abundance','disease-associated host state','malignant epithelial state','sampling / treatment / technical confounding'],'status':'active'})
urls={
 'ascites2023':'https://www.nature.com/articles/s41467-023-36310-9',
 'spatial2024':'https://www.gastrojournal.org/article/S0016-5085(24)05348-4/fulltext',
 'cgnl_ko2023':'https://www.mdpi.com/2073-4409/12/15/1981',
 'cgnl_ko_europepmc':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10417749/fullTextXML',
 'ascites_europepmc':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9929081/fullTextXML',
 'GSE308231':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308231&targ=self&form=text&view=full',
 'GSE308231_files':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308231/suppl/',
 'spatial_search':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22Spatially%20Resolved%20Niche%22&format=json',
}
def fetch(k,u):
 p=N/'literature'/(k+'.response');r={'key':k,'url':u}
 try:
  a=requests.get(u,timeout=(15,45));p.write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest(),path=str(p),content_type=a.headers.get('Content-Type'));print(k,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=str(e);print(k,str(e),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:receipts=list(ex.map(lambda kv:fetch(*kv),urls.items()))
put(N/'metadata/ACQUISITION_01.json',receipts)
print('ROOT',str(N))
