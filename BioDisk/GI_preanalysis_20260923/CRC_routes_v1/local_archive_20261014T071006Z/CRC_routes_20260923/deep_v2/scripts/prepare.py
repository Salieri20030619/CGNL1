import pathlib,json,csv,urllib.request,concurrent.futures,time
B=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');O=B/'deep_v2';rows=list(csv.DictReader((B/'results/route1/cross_subtype_sensitivity_screen.tsv').open(),delimiter='\t'));genes=[r['gene'] for r in rows if r['passes_descriptive_sensitivity']=='True']
p={'date':'2026-09-24','frozen_genes':genes,'goal':'A credible mechanism-oriented project proposal for approximately IF10 target, conditional on independent bioinformatics evidence and subsequent wet validation, no publication guarantee.','analyses':['Patient-specific RNA-inferred CNV support with nonmalignant fibroblast/endothelial references; exclude all126 candidate genes from CNV inference; predefined95th/99th reference quantile sensitivities; do not equate RNA-inferred CNV with DNA confirmation','Independent GSE50760 paired primary-liver validation: genomewide main FDR and fixed126 validation family reported separately, with liver/stromal/epithelial composition sensitivities; site eligibility audited','Cross-source localization and spatial support for data-supported candidates; not require all mechanisms to peak at primary invasive front','Focused novelty check of candidates meeting data gates; avoid common already-studied liver-like reprogramming axes','Specify discriminating wet experiments and stop/go criteria for top1-2 mechanisms rather than generic proliferation/EMT experiments'],'statistics':'Patient biological unit; discovery five pairs remains exploratory. Independent bulk n>=10 reports effects,P,FDR. Families separated by genomewide DE, frozen candidate validation, program correlation and sensitivity. No significance rescue by dropping failing genes.'}
(O/'PLAN.json').write_text(json.dumps(p,indent=2))
jobs=[('GSE50760.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/soft/GSE50760_family.soft.gz'),('GSE50760.files.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/suppl/'),('refseq_hg38.txt.gz','https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz'),('infercnvpy_pypi.json','https://pypi.org/pypi/infercnvpy/json')]
def fetch(v):
 name,url=v;op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}));t=time.time()
 try:
  with op.open(url,timeout=35) as r:data=r.read()
  (O/'metadata'/name).write_bytes(data);j={'file':name,'url':url,'bytes':len(data),'exit_code':0}
 except Exception as e:j={'file':name,'url':url,'exit_code':1,'error':repr(e)}
 print(json.dumps(j),flush=True);return j
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(fetch,jobs))
(O/'metadata_fetch.json').write_text(json.dumps(r,indent=2))
