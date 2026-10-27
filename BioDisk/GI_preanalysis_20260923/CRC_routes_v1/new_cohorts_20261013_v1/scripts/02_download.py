from pathlib import Path
import requests,concurrent.futures,json,hashlib,gzip,datetime,re
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'new_cohorts_20261013_v1'
for p in [R/'scripts',R/'metadata',R/'raw',R/'logs',R/'results',R/'figures']:p.mkdir(parents=True,exist_ok=True)
old=json.loads((B/'increment_context_20261012_v1/SPECIFICATION.json').read_text())
spec=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),target='CGNL1',new_cohort_target_results_seen=False,
primary_cohorts=['GSE200997','GSE132465'],primary_cohort_role='new primary-CRC epithelial information validation; not liver colonization validation',
liver_descriptive=['GSE178318','GSE231559'],liver_role='Previously downloaded datasets, new epithelial endpoint; small cohorts, eligibility and descriptive results retained; cannot count reanalyses as new patients',
panels=old['panels'],outcomes={'actomyosin':old['outcomes']['actomyosin'],'cycling':old['panels']['cycling']},
ordinary_markers=old['ordinary_markers'],
selection='Tumor samples only; author epithelial annotation as primary, frozen marker intersection sensitivity; retain CGNL1-zero cells; >=20 epithelial cells per patient; fixed QC UMI>=500, genes>=200, mito<=25%; gene/program coverage>=80%; all7 actomyosin genes required',
aggregation='sum raw counts per patient, ln(1+CP10k) using all measured genes as denominator; record >=20 and >=50-cell coverage, latter only prespecified sensitivity',
models={'M0':['depth'],'M1':['depth','epithelial_identity','mature_colonocyte'],'M1_VIL1':['depth','epithelial_identity','mature_colonocyte','VIL1'],'M1_KRT20':['depth','epithelial_identity','mature_colonocyte','KRT20'],'M1_all6':['depth','epithelial_identity','mature_colonocyte']+old['ordinary_markers']},
inference='Patient biological unit. HC3 t,full rank,residual_df>=5. n<10 descriptive coefficient/spearman only in this frozen validation; no inferred disease association from cells. Raw targetP<.05 nominal hard gate; BH across2 endpoints within dataset/model/definition; keep FDR and exact family. No screening by FDR alone.',
prediction='LOPO baseline vs addedCGNL1; all preprocessing training only; no tuning. Cross-cohort transfer original23->each new cohort in allgene units; distribution mismatch disclosed. Ordinary six comparator LOPO benchmark with same M1.',
decision='Do both original actomyosin and exploratory cycling provide consistent independent new-patient increments? Separate primary CRC generalization from liver-specific evidence. No target switch or changing signatures from result.',
lineage_markers=json.loads((B/'formal_v1_20261002/frozen/EXECUTION_CONFIG.json').read_text())['lineage_markers'])
p=R/'SPECIFICATION.json'
if not p.exists():
 p.write_text(json.dumps(spec,indent=2));(R/'SPEC_FREEZE.json').write_text(json.dumps(dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),utc=spec['created_utc']),indent=2))
jobs=[]
for ds,ext in [('GSE200997','csv'),('GSE132465','txt')]:
 stem='https://ftp.ncbi.nlm.nih.gov/geo/series/'+ds[:-3]+'nnn/'+ds+'/suppl/'
 for suffix in ['cell_annotation','raw_UMI_count_matrix']:
  name=ds+'_GEO_processed_CRC_10X_'+suffix+'.'+ext+'.gz';jobs.append((name,stem+name))
 jobs.append((ds+'_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/'+ds[:-3]+'nnn/'+ds+'/soft/'+ds+'_family.soft.gz'))
jobs.append(('atlas_tree.json','https://api.github.com/repos/icbi-lab/crc-atlas/git/trees/main?recursive=1'))
def fetch(item):
 name,url=item;dest=(R/'raw' if 'matrix' in name else R/'metadata')/name
 rec=dict(name=name,url=url,path=str(dest),attempts=[])
 for attempt in range(2):
  try:
   if not dest.exists():
    sess=requests.Session();sess.trust_env=False
    with sess.get(url,stream=True,timeout=(20,90)) as response:
     response.raise_for_status();size=response.headers.get('Content-Length');h=hashlib.sha256();n=0;tmp=dest.with_suffix(dest.suffix+'.part')
     with tmp.open('wb') as f:
      for block in response.iter_content(1048576):f.write(block);h.update(block);n+=len(block)
     if size:assert n==int(size)
     tmp.replace(dest)
   if name.endswith('.gz'):
    with gzip.open(dest,'rb') as f:
     while f.read(1048576):pass
   rec.update(status='verified',bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
   break
  except Exception as e:rec['attempts'].append(repr(e));rec['status']='failed'
 print('FETCH',name,rec['status'],rec.get('bytes'),flush=True)
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(fetch,jobs))
(R/'DOWNLOAD_MANIFEST.json').write_text(json.dumps(records,indent=2))
(R/'DOWNLOAD_COMPLETE.json').write_text(json.dumps(dict(all_success=all(x['status']=='verified' for x in records),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
for ds in ['GSE200997','GSE132465']:
 p=next((R/'metadata').glob(ds+'*annotation*'),None)
 if p:
  with gzip.open(p,'rt') as f:
   for i in range(4):print('ANNOTATION',ds,f.readline().strip()[:1500],flush=True)
print('DOWNLOAD_PHASE_COMPLETE',flush=True)
