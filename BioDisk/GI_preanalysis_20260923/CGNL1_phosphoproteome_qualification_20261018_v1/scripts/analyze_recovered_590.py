from pathlib import Path
import csv,json,re,math,collections,statistics,hashlib,datetime
import numpy as np
import openpyxl
from scipy.stats import rankdata,wilcoxon
csv.field_size_limit(10000000)
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576')
FAMILY=['S112','S198','S202','S261','S283','S296','S297','S349','S414','S486','S514','S708']
SEQ=''.join(x.strip() for x in (B/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not x.startswith('>'))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def signed(x):
 y=[v for v in x if v!=0]
 if not y:return None
 ranks=[int(round(2*v)) for v in rankdata(np.abs(y))];obs=sum(r for r,v in zip(ranks,y) if v>0);cut=min(obs,sum(ranks)-obs)
 ways=[0]*(cut+1);ways[0]=1
 for r in ranks:
  for k in range(cut,r-1,-1):ways[k]+=ways[k-r]
 return min(1,2*sum(ways)/2**len(y))
def load(path,mode,phos):
 with path.open(encoding='utf-8') as f:
  reader=csv.DictReader(f,delimiter='\t');heads=reader.fieldnames
  cols=[c for c in heads if (c.startswith('Log ') if mode=='B03' else re.fullmatch('Intensity [0-9]+_?[PT]',c))]
  arrays=[];targets=[]
  for row in reader:
   if row[heads[0]] in ['Mean','Median','StdDev'] or row.get('Reverse')=='+' or row.get('Potential contaminant')=='+':continue
   a=[]
   for c in cols:
    try:v=float(row[c])
    except (ValueError,TypeError):v=np.nan
    if not math.isfinite(v) or (mode=='H02' and v<=0):v=np.nan
    elif mode=='H02':v=math.log2(v)
    a.append(v)
   arrays.append(a)
   if row.get('Gene')=='CGNL1' or 'CGNL1' in row.get('Gene names','').split(';'):targets.append((len(arrays)-1,row))
 matrix=np.asarray(arrays);offset=np.nanmedian(matrix,axis=0)
 # Verify independent standard-library medians for every sample column.
 for i in range(len(cols)):
  v=matrix[:,i];v=v[np.isfinite(v)]
  if len(v):assert abs(statistics.median(v.tolist())-offset[i])<1e-12
 labels=[c[4:].split(':')[0] if mode=='B03' else c[10:].replace('_','') for c in cols]
 lookup=collections.defaultdict(list)
 for i,label in enumerate(labels):lookup[label].append(i)
 return matrix,offset,targets,lookup,cols
def run(cohort):
 o=B/('results/'+cohort+'_recovered_590');o.mkdir(exist_ok=True)
 spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),family=FAMILY,primary='Within clinically mapped donor: median-centered phosphosite T-N difference minus median-centered independent total-protein T-N difference',sensitivity='Same paired endpoint without extra sample-median centering',small_n='n<10 descriptive per frozen project specification',missing='No new imputation; raw zero/nonpositive intensity excluded',normalization='Each omics separately; all eligible source rows; H02 raw intensity log2, B03 source Log columns retained',localization='H02 global localization probability >=0.75, unique sequence-window map; B03 singly modified source peptides uniquely mapped to canonical sequence',inference='Two-sided exact signed rank, BH among evaluable prespecified sites within this cohort',limitations=['B03 source log base not independently confirmed','H02 MaxQuant match-between-runs enabled; do not describe every positive intensity as separately localized MS/MS evidence'])
 save(o/'SPEC.json',spec)
 if cohort=='H02':
  protein=next((D/cohort/'global_589').rglob('proteinGroups.txt'));phospho=D/cohort/'selected_585/Phospho (STY)Sites.txt'
  w=openpyxl.load_workbook(next((D/cohort/'global_589').rglob('*.xlsx')),read_only=True,data_only=True);mapping=collections.defaultdict(set)
  for sheet in w:
   if not sheet.title.startswith('Quantitative analysis'):continue
   for row in sheet.iter_rows(min_row=2,values_only=True):
    if len(row)>=5 and row[2] and row[4] and str(row[3]) in ['A','0']:mapping[str(row[2]).replace('_','')].add(str(row[4]))
  assert all(len(v)==1 for v in mapping.values());mapping={k:next(iter(v)) for k,v in mapping.items()}
 else:
  protein=D/cohort/'iCC_NCC_Proteome.tmt10.tsv';phospho=D/cohort/'iCC_NCC_Phosphoproteome.phosphosite.tmt10.tsv'
 pa,po,pt,pl,pc=load(protein,cohort,False);sa,so,st,sl,sc=load(phospho,cohort,True);assert len(pt)==1
 pairs=[];excluded=[]
 if cohort=='H02':
  rawpairs=sorted(k[:-1] for k in sl if k.endswith('T') and k[:-1]+'P' in sl)
  for case in rawpairs:
   t,n=case+'T',case+'P'
   if t not in mapping or n not in mapping or t not in pl or n not in pl:excluded.append(dict(raw=case,reason='Outside source early-stage BCLC A/0 quantitative cohort, missing clinical crosswalk, or missing independent total-protein column'));continue
   assert mapping[t].endswith('_T') and mapping[n].endswith('_N') and mapping[t][:-2]==mapping[n][:-2]
   pairs.append(dict(case=mapping[t][:-2],tumor=t,normal=n,raw=case))
 else:
  clinical=[]
  for sid in ['PDC000356','PDC000363']:
   rows=json.loads((B/('metadata/bulk_576/'+sid+'_biospecimens.json')).read_text())['biospecimenPerStudy'];cases=collections.defaultdict(dict)
   for r in rows:
    if r['sample_type']=='Normal Adjacent Tissue':role='normal'
    elif r['sample_type']=='Human Tumor Original Cells':role='tumor'
    else:continue
    assert role not in cases[r['case_submitter_id']]
    cases[r['case_submitter_id']][role]=r['sample_submitter_id']
   clinical.append(dict(cases))
  assert clinical[0]==clinical[1]
  for case,v in clinical[0].items():
   if len(v)==2 and all(v[k] in pl and v[k] in sl for k in ['tumor','normal']):pairs.append(dict(case=case,**v))
   else:excluded.append(dict(case=case,reason='Incomplete paired matrix coverage'))
 assert len({p['case'] for p in pairs})==len(pairs)
 save(o/'PAIRING.json',dict(eligible=pairs,excluded=excluded))
 mapped=collections.defaultdict(list);mappingqa=[]
 for i,row in st:
  positions=[]
  if cohort=='H02':
   if float(row['Localization prob'])<0.75:continue
   win=row['Sequence window'].split(';')[0];unpad=win.strip('_');starts=[m.start() for m in re.finditer('(?='+re.escape(unpad)+')',SEQ)]
   if len(starts)==1:positions=[row['Amino acid']+str(starts[0]+16-(len(win)-len(win.lstrip('_'))))]
  else:
   for peptide in row['Peptide'].split(';'):
    starts=[m.start() for m in re.finditer('(?='+re.escape(peptide.upper())+')',SEQ)]
    if len(starts)!=1:continue
    sites=[aa.upper()+str(starts[0]+j+1) for j,aa in enumerate(peptide) if aa.islower()]
    if len(sites)==1:positions.extend(sites)
  sites=set(positions);mappingqa.append(dict(source=row.get('Phosphosite',row.get('Position')),mapped=list(sites)))
  if len(sites)==1:mapped[next(iter(sites))].append(i)
 save(o/'SEQUENCE_QA.json',mappingqa)
 save(o/'NORMALIZATION.json',dict(protein=dict(zip(pc,po.tolist())),phosphosite=dict(zip(sc,so.tolist()))))
 pi=pt[0][0];values=[];results=[]
 def sample(matrix,offset,row,label,lookup,center):
  a=[matrix[row,k]-(offset[k] if center else 0) for k in lookup[label] if np.isfinite(matrix[row,k])]
  return float(np.median(a)) if a else None
 for site in FAMILY:
  effect=[];raw=[]
  for pair in pairs:
   endpoints=[]
   for centered in [True,False]:
    pvals=[sample(pa,po,pi,pair[k],pl,centered) for k in ['tumor','normal']]
    ds=[]
    for j in mapped.get(site,[]):
     v=[sample(sa,so,j,pair[k],sl,centered) for k in ['tumor','normal']]
     if all(x is not None for x in v):ds.append(v[0]-v[1])
    endpoints.append(float(np.median(ds))-pvals[0]+pvals[1] if ds and all(x is not None for x in pvals) else None)
   if endpoints[0] is None:continue
   effect.append(endpoints[0]);raw.append(endpoints[1]);values.append(dict(site=site,case=pair['case'],adjusted=endpoints[0],uncentered=endpoints[1]))
  n=len(effect);p=signed(effect) if n>=10 else None;pr=signed(raw) if n>=10 else None
  if p is not None and len(set(abs(v) for v in effect))==n and all(v!=0 for v in effect):assert abs(p-float(wilcoxon(effect,method='exact').pvalue))<1e-12
  results.append(dict(site=site,n=n,up=sum(v>0 for v in effect),down=sum(v<0 for v in effect),mean=statistics.mean(effect) if n else None,median=statistics.median(effect) if n else None,P=p,uncentered_mean=statistics.mean(raw) if n else None,uncentered_P=pr))
 tested=sorted([r for r in results if r['P'] is not None],key=lambda r:r['P']);q=1
 for i in range(len(tested)-1,-1,-1):q=min(q,tested[i]['P']*len(tested)/(i+1));tested[i]['FDR']=q
 save(o/'PATIENT_VALUES.json',values);save(o/'RESULT.json',dict(cohort=cohort,eligible_donors=len(pairs),results=results,family_size=len(tested),status='exploratory_complete_pending_independent_review'))
 save(o/'SOURCE_MANIFEST.json',[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [protein,phospho]])
 print(cohort,'eligible',len(pairs),json.dumps(results),flush=True)
for cohort in ['B03','H02']:run(cohort)
