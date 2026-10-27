from pathlib import Path
import json,gzip,csv,math,statistics,hashlib,datetime,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C04_full_fraction_qualification/reporter_normalization_20261026')
channels=['127n','127c','128n','128c','129n','129c','130n','130c']
def save(n,x):(P/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(n,fields,rows):
 with (P/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
assert json.loads((P/'EXTRACTION_ACCEPTANCE.json').read_text())['pass_all']
assert not (P/'NORMALIZATION_FACTORS.json').exists()
save('FORMULA_SPEC_BEFORE_FACTORS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),channel_order=channels,aggregation='sum reference-eligible unique MS3 intensities within plain peptide plus phosphate count across fractions',N0='log2 offset 0; raw diagnostic',N1='For each feature log2(I) minus mean across8channels, then channel median over features, then subtract mean of8channel medians; corrected_log2=log2(I)-offset',N2='log2(channel summed reference feature intensity) minus mean of8channel log2 totals; corrected_log2=log2(I)-offset',Cgnl1_excluded=True,missing_not_imputed=True,source_exact_reproduction=False,biological_tests=0))
receipts=json.loads((P/'EXTRACTION_EXIT.json').read_text())['files'];assert len(receipts)==12
features={};reference=[];targets=[];seen=set()
for receipt in receipts:
 path=P/(Path(receipt['file']).stem+'.reporters.json.gz');assert hashlib.sha256(path.read_bytes()).hexdigest()==receipt['output_sha256']
 with gzip.open(path,'rt',encoding='utf-8') as f:records=json.load(f)['records']
 for x in records:
  key=(x['file'],x['MS3']);assert key not in seen;seen.add(key)
  if x['reference_eligible']:
   assert not x['Cgnl1_associated'] and x['all8_positive_unique'] and not x['assignment_ambiguous']
   feat=tuple(x['feature']);v=[r['intensity'] for r in x['reporters']];assert len(v)==8 and all(a>0 for a in v)
   if feat not in features:features[feat]=[0.0]*8
   features[feat]=[a+b for a,b in zip(features[feat],v)]
   reference.append(dict(file=x['file'],MS3=x['MS3'],peptide=feat[0],phosphate_count=feat[1],**dict(zip(channels,v))))
  if x['Cgnl1_associated']:targets.append(x)
assert len(reference)==5068 and len(targets)==4
feature_rows=[dict(peptide=k[0],phosphate_count=k[1],**dict(zip(channels,v))) for k,v in sorted(features.items())]
table('REFERENCE_MS3.tsv',['file','MS3','peptide','phosphate_count']+channels,reference)
table('REFERENCE_FEATURES.tsv',['peptide','phosphate_count']+channels,feature_rows)
centered=[]
for v in features.values():
 log=[math.log2(z) for z in v];m=statistics.mean(log);centered.append([z-m for z in log])
n1=[statistics.median(row[i] for row in centered) for i in range(8)];mean=statistics.mean(n1);n1=[v-mean for v in n1]
n2=[math.log2(sum(v[i] for v in features.values())) for i in range(8)];mean=statistics.mean(n2);n2=[v-mean for v in n2]
offsets={'N0':[0.0]*8,'N1':n1,'N2':n2};factors=[dict(model=m,channel=c,log2_offset=v,divisor=2**v) for m,vs in offsets.items() for c,v in zip(channels,vs)]
save('NORMALIZATION_FACTORS.json',dict(factors=factors,reference_MS3=len(reference),reference_features=len(features),Cgnl1_excluded=True,whole_project_complete=False))
table('NORMALIZATION_FACTORS.tsv',['model','channel','log2_offset','divisor'],factors)
raw=[];normalized=[]
for x in targets:
 for i,c in enumerate(channels):
  v=x['reporters'][i]['intensity'];raw.append(dict(file=x['file'],MS3=x['MS3'],channel=c,intensity=v,matches=x['reporters'][i]['matches']))
  for model,off in offsets.items():normalized.append(dict(file=x['file'],MS3=x['MS3'],channel=c,model=model,intensity=v,normalized_log2=math.log2(v)-off[i] if v is not None else None))
table('CGNL1_RAW_REPORTERS.tsv',['file','MS3','channel','intensity','matches'],raw)
table('CGNL1_NORMALIZED_REPORTERS.tsv',['file','MS3','channel','model','intensity','normalized_log2'],normalized)
save('NORMALIZATION_RESULT.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reference_MS3=len(reference),reference_features=len(features),target_MS3=len(targets),target_observed_channels=sum(x['intensity'] is not None for x in raw),target_missing_or_ambiguous_channels=sum(x['intensity'] is None for x in raw),factor_count=len(factors),model_N1_N2_max_log2_offset_difference=max(abs(a-b) for a,b in zip(n1,n2)),independent_QA_pending=True,disease_tests=0,protein_adjustment_not_performed=True,whole_project_complete=False))
print(json.dumps(dict(reference_MS3=len(reference),reference_features=len(features),factor_count=len(factors),N1=n1,N2=n2)))
