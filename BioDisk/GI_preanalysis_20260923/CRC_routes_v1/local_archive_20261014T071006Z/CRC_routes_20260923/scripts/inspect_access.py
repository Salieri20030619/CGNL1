import pathlib,json,gzip,csv,re,subprocess,concurrent.futures,sys
sys.stdout.reconfigure(encoding='utf-8');O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');M=O/'metadata';j=json.loads((M/'zenodo_latest.json').read_text());print('LATEST',j.get('id'),j.get('metadata',{}).get('access_right'),[(x.get('key'),x.get('size')) for x in j.get('files',[])])
j=json.loads((M/'github_tree.json').read_text());print('GITHUB',[(x['path'],x.get('size')) for x in j.get('tree',[])])
s=(M/'spatial2026.html').read_text(encoding='utf-8');links=re.findall(r'href="([^"]+)"',s);print('SUPP',[x for x in links if 'supplement' in x or '.xlsx' in x][:20]);p=json.loads((M/'GSE225857_paper.xml').read_text());print('PAPER',p['resultList']['result'][0].get('pmcid'))
with gzip.open(M/'nonimmune_meta.txt.gz','rt') as f:
 for _ in range(4):print(f.readline()[:1800])
pmcid=p['resultList']['result'][0]['pmcid'];jobs=[('nonimmune_counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7058nnn/GSM7058755/suppl/GSM7058755_non_immune_counts.txt.gz'),('GSE225857_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML')]
# Start article and small supplementary downloads first; large counts run separately.
jobs=jobs[1:]+[(x.split('/')[-1],('https://pmc.ncbi.nlm.nih.gov'+x if x.startswith('/') else 'https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/'+x)) for x in links if x.endswith('.xlsx')][:2]
def run(v):
 n,u=v;r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','120','-o',str(M/n),u],capture_output=True);return n,r.returncode,(M/n).stat().st_size if (M/n).exists() else 0
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:print(list(ex.map(run,jobs)))
