import pathlib,json,gzip,subprocess,concurrent.futures
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');j=json.loads((O/'metadata/zenodo19043057.json').read_text());print(json.dumps({k:v for k,v in j.items() if k!='links'},ensure_ascii=False)[:6500])
s=gzip.open(O/'metadata/GSE231559.soft.gz','rt',encoding='utf-8').read()
for c in s.split('^SAMPLE = ')[1:]:print(c.splitlines()[0],'; '.join(x for x in c.splitlines() if x.startswith(('!Sample_title','!Sample_characteristics_ch1'))))
jobs=[('zenodo_files.json',j['links']['files']),('zenodo_latest.json',j['links']['latest']),('github_tree.json','https://api.github.com/repos/yliuup/CRC_micromets_ST/git/trees/main?recursive=1'),('spatial2026.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/'),('GSE225857_paper.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:37327339&format=json'),('nonimmune_meta.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7058nnn/GSM7058755/suppl/GSM7058755_non_immune_meta.txt.gz')]
def run(v):
 n,u=v;r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','90','-o',str(O/'metadata'/n),u],capture_output=True);return n,r.returncode,(O/'metadata'/n).stat().st_size if (O/'metadata'/n).exists() else 0
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:print(list(e.map(run,jobs)))
