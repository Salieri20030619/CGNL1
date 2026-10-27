from pathlib import Path
import csv,json,re,math,statistics,collections,datetime
from scipy.stats import wilcoxon
csv.field_size_limit(10000000)
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576')
seq=''.join(s.strip() for s in (B/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not s.startswith('>'))
for cohort in ['B03','H02']:
 o=B/('results/'+cohort+'_recovered_590');res=json.loads((o/'RESULT.json').read_text());pairs=json.loads((o/'PAIRING.json').read_text())['eligible'];offset=json.loads((o/'NORMALIZATION.json').read_text());values=json.loads((o/'PATIENT_VALUES.json').read_text())
 manifest=json.loads((o/'SOURCE_MANIFEST.json').read_text());targets=[];lookups=[]
 for source in manifest:
  with Path(source['path']).open() as f:
   reader=csv.DictReader(f,delimiter='\t');lookup={}
   for col in reader.fieldnames:
    if cohort=='B03' and col.startswith('Log '):label=col[4:].split(':')[0]
    elif cohort=='H02' and re.fullmatch('Intensity [0-9]+_?[PT]',col):label=col[10:].replace('_','')
    else:continue
    if label=='ref01':continue
    assert label not in lookup;lookup[label]=col
   rows=[r for r in reader if r.get('Gene')=='CGNL1' or 'CGNL1' in r.get('Gene names','').split(';')]
  targets.append(rows);lookups.append(lookup)
 assert len(targets[0])==1
 site_rows=collections.defaultdict(list)
 for row in targets[1]:
  sites=set()
  if cohort=='H02':
   if float(row['Localization prob'])<.75:continue
   win=row['Sequence window'].split(';')[0];pattern=win.strip('_');start=seq.find(pattern)
   if start>=0 and seq.find(pattern,start+1)<0:sites.add(row['Amino acid']+str(start+16-(len(win)-len(win.lstrip('_')))))
  else:
   for peptide in row['Peptide'].split(';'):
    marks=[i for i,c in enumerate(peptide) if c.islower()];start=seq.find(peptide.upper())
    if len(marks)==1 and start>=0 and seq.find(peptide.upper(),start+1)<0:sites.add(peptide[marks[0]].upper()+str(start+marks[0]+1))
  if len(sites)==1:site_rows[next(iter(sites))].append(row)
 expected={(v['site'],v['case']):v for v in values};errors=[];maxerr=0;allvalues=collections.defaultdict(list)
 def difference(row,pair,lookup,off,center):
  t,n=lookup[pair['tumor']],lookup[pair['normal']]
  try:tv,nv=float(row[t]),float(row[n])
  except (ValueError,TypeError):return None
  if not math.isfinite(tv) or not math.isfinite(nv):return None
  if cohort=='H02':
   if tv<=0 or nv<=0:return None
   value=math.log2(tv)-math.log2(nv)
  else:value=tv-nv
  return value-off[t]+off[n] if center else value
 for endpoint in res['results']:
  for pair in pairs:
   computed=[]
   for center in [True,False]:
    dp=difference(targets[0][0],pair,lookups[0],offset['protein'],center)
    ds=[difference(row,pair,lookups[1],offset['phosphosite'],center) for row in site_rows.get(endpoint['site'],[])];ds=[v for v in ds if v is not None]
    computed.append(statistics.median(ds)-dp if ds and dp is not None else None)
   if computed[0] is None:
    assert (endpoint['site'],pair['case']) not in expected;continue
   ref=expected[(endpoint['site'],pair['case'])];maxerr=max(maxerr,abs(computed[0]-ref['adjusted']),abs(computed[1]-ref['uncentered']));allvalues[endpoint['site']].append(computed[0])
  v=allvalues[endpoint['site']];assert len(v)==endpoint['n']
  if len(v)>=10:
   assert len(set(abs(z) for z in v))==len(v) and all(z!=0 for z in v),'Exact scipy crosscheck requires untied nonzero values'
   p=float(wilcoxon(v,method='exact').pvalue);assert abs(p-endpoint['P'])<1e-12
  else:assert endpoint['P'] is None
 tested=[x for x in res['results'] if x['P'] is not None]
 for row in tested:
  # Independent BH definition: min over all larger sorted p values.
  ordered=sorted(x['P'] for x in tested);idx=ordered.index(row['P']);q=min([1]+[len(ordered)*ordered[k]/(k+1) for k in range(idx,len(ordered))]);assert abs(q-row['FDR'])<1e-12
 assert maxerr<1e-12
 qa=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_pass=True,cohort=cohort,patient_site_rows=len(values),max_effect_error=maxerr,scope='Separate direct source-target reconstruction, sequence mapping, eligibility counts, scipy exact statistics and BH verification. Normalization offsets reused; original calculation also checked every offset against stdlib medians. Clinical interpretation and source design not proved by numerical QA.')
 (o/'INDEPENDENT_QA591.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa),flush=True)
