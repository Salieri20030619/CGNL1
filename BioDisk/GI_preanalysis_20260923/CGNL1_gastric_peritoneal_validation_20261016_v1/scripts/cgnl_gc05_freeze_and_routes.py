from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,subprocess,re
from concurrent.futures import ThreadPoolExecutor
from bs4 import BeautifulSoup
O=Path.cwd();B=O.parent
spec={'utc':datetime.now(timezone.utc).isoformat(),'target':'CGNL1','target_results_not_yet_examined_in_bulk_or_sc':True,'primary_bulk':'GSE314812 17 previously qualified matched patients, GCPM minus GCPT in existing independently documented TMM logCPM; author RSEM expected counts not UMI; validate input against author file','bulk_estimand':'site-associated tissue difference, distinct from2022 primary-tumor P+/Cy+ vs P0/Cy0 contrast','bulk_models':{'M0':'paired change mean and exact Wilcoxon; main single target; also t/HC3 estimate','M1':'paired change ~ intercept + paired FFPE-indicator difference; target intercept at matched preservation','M2':'M1 + paired changes in epithelial, stroma, myeloid, adipocyte, endothelial, mesothelial proxies; intercept is conditional site difference; full rank and residual df>=5 required'},'bulk_FDR':'three HC3 target-intercept tests same sensitivity family; primary Wilcoxon single prespecified endpoint no multiple-gene FDR; all estimates/P preserved','bulk_robustness':['leave-one-patient mean range','same preservation pairs','naive full-library logCPM sensitivity from author expected counts;no arbitrary normalizationselection'],'sc_source':'GSE228598 target-excluded frozen lineage gates reused, >=20 cells/patient/lineage, preserve >=50 sensitivity; raw counts sum and detection per lineage; no ascites-versus-washings inferential claim withoutclinicalqualification','sc_paired':'GSE234129 author subtypes and patient mapping; no forced malignant relabel; source composition and within-lineage Ca/As differences descriptive if<10 independent pairs','normalization':'singlecell patient-lineage pooled ln1pCP10k using all-gene raw UMI; preserve zero cells and no unmeasured-as-zero','GSE239676_gate':'downloaded feature list8630 rows contains neitherCGNL1 norENSG00000128849; metadata download partial and not analysis input; do not downloadlargefilteredmatrixasCGNL1validation','original2022':'publicfilelisting not yet resolved; continue official repository diagnosis, no reconstruction of patient expression fromgroupfigure','status':'frozen_before_target_results'}
b=(json.dumps(spec,ensure_ascii=False,indent=2)+'\n').encode('utf-8');(O/'ANALYSIS_SPECIFICATION.json').write_bytes(b);(O/'SPEC_FREEZE.json').write_text(json.dumps({'utc':spec['utc'],'sha256':hashlib.sha256(b).hexdigest()},indent=2),encoding='utf-8')
# Mark earlier incomplete metadata inspections as provisional, not valid cohort counts.
(O/'metadata/GSE239676_METADATA_STATUS.json').write_text(json.dumps({'status':'incomplete_download_not_used_for_analysis','earlier_printed_category_counts':'partial-file diagnostics only, not qualified sample sizes','feature_file':'transport complete, plain text despite .gz suffix; target absent; no zero-expression claim'},indent=2),encoding='utf-8')
soup=BeautifulSoup((O/'metadata/mendeley_page').read_text(encoding='utf-8'),'html.parser');scripts=[s.get('src') for s in soup.find_all('script') if s.get('src')];print('JS',scripts)
urls={}
for i,u in enumerate(scripts):
 if '/datasets/' in u and '.js' in u:urls['mendeley_js_'+str(i)]=u if u.startswith('http') else 'https://data.mendeley.com'+u
urls['GSE314812_author_download']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314812/suppl/GSE314812_rsem.merged.gene_counts.tsv.gz'
def fetch(kv):
 name,u=kv;p=(O/'raw/GSE314812_author_counts.download') if 'author_download' in name else O/'metadata'/name
 r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','15','--max-time','600' if 'author_download' in name else '45','--retry','1','-sS','-w','%{http_code}','-o',str(p),u],capture_output=True,timeout=1230 if 'author_download' in name else 110);rec={'name':name,'url':u,'exit':r.returncode,'http':r.stdout.decode('utf-8'),'stderr':r.stderr.decode('utf-8')}
 if p.exists():
  raw=p.read_bytes();rec.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
  if '.js' in u:
   txt=raw.decode('utf-8');rec['route_excerpts']=[txt[max(0,m.start()-220):m.start()+350] for m in list(re.finditer(r'folder_id|/files|/folders',txt))[:20]]
 print(json.dumps(rec),flush=True);return rec
with ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(fetch,urls.items()))
(O/'metadata/DOWNLOAD_RECEIPTS_05.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
# concise phenotype/schema evidence
import pandas as pd
meta=pd.read_csv(B/'02_raw/GSE234129_meta.tsv.gz',sep='\t');print('234_CATS',{c:meta[c].value_counts().to_dict() for c in meta if meta[c].nunique()<40})
m=pd.read_csv(B/'GC_validation_v2/candidate_discovery/anatomy_proxy_sample_scores.tsv',sep='\t');print('PRESERVATION',pd.crosstab(m.site,m.preservation).to_dict());print('BULK_COLUMNS',list(m))
