from pathlib import Path
import json,csv,hashlib,datetime,math,statistics
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');F=O/'C04_full_fraction_qualification';N=F/'reporter_normalization_20261026';W=F/'disease_descriptive_20261026';W.mkdir(exist_ok=True)
def save(n,x):(W/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert json.loads((N/'NORMALIZATION_ACCEPTANCE.json').read_text())['independent_QA_pass']
assert not (W/'SPEC_BEFORE_GROUP_COMPARISON.json').exists()
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),question='Inflamed versus control relative Cgnl1 phosphopeptide abundance, before total-protein adjustment',designation='exploratory reanalysis; earlier pilot abundance data and study groups already known, not blinded preregistration',biological_n=8,groups={'control':5,'inflamed':3},paired=False,endpoints=[dict(id='D01',role='primary measurement-qualified candidate',file='phospho_fraction_11.raw',scans=[19264],identity='mouse S284; human S283 sequence correspondence; not human disease measurement'),dict(id='D02',role='secondary unresolved-site phosphopeptide abundance',file='phospho_fraction_1.raw',scans=[26697,26749],identity='SPTAPSPQAYSETK with one phosphate; site-isomer ambiguity retained; sum unique complete MS3 intensities per animal before log transform'),dict(id='D03',role='incomplete S284 measurement diagnostic',file='phospho_fraction_11.raw',scans=[19300],identity='three multiple-centroid channels unavailable; preserve observed values')],models={'N0':'uncorrected diagnostic','N1':'primary exploratory median-feature normalization','N2':'summed-reference-intensity sensitivity'},estimands=['mean log2 inflamed minus mean log2 control','difference of group medians on log2 scale','individual inflamed values versus control median','leave-one-animal-out mean difference range'],minimum_evaluable_per_group=2,missing_values='no imputation; mark comparison not evaluable if minimum unmet',inference='descriptive only, no P, FDR, CI or biological replicate inflation',no_model_or_endpoint_selected_by_direction=True,total_protein_adjustment=False,causality_or_occupancy=False,inputs={n:sha(N/n) for n in ['CGNL1_RAW_REPORTERS.tsv','NORMALIZATION_FACTORS.tsv','NORMALIZATION_ACCEPTANCE.json']},mapping_sha256=sha(O/'C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json'))
save('SPEC_BEFORE_GROUP_COMPARISON.json',spec)
mapping=json.loads((O/'C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json').read_text())['records'];assert len(mapping)==8 and len({x['sample'] for x in mapping})==8
assert {g:sum(x['group']==g for x in mapping) for g in ['control','inflamed']}==spec['groups']
with (N/'CGNL1_RAW_REPORTERS.tsv').open() as f:raw={(x['file'],int(x['MS3']),x['channel']):float(x['intensity']) if x['intensity'] else None for x in csv.DictReader(f,delimiter='\t')}
with (N/'NORMALIZATION_FACTORS.tsv').open() as f:offsets={(x['model'],x['channel']):float(x['log2_offset']) for x in csv.DictReader(f,delimiter='\t')}
values=[];results=[];loo=[]
for endpoint in spec['endpoints']:
 for model in spec['models']:
  rows=[]
  for animal in mapping:
   c=animal['phospho_TMT_channel'];obs=[raw[(endpoint['file'],scan,c)] for scan in endpoint['scans']];abundance=sum(obs) if all(v is not None for v in obs) else None
   value=math.log2(abundance)-offsets[(model,c)] if abundance is not None else None
   row=dict(endpoint=endpoint['id'],model=model,sample=animal['sample'],group=animal['group'],channel=c,source_MS3_count=len(obs),raw_summed_intensity=abundance,log2_value=value);values.append(row);rows.append(row)
  ctrl=[x['log2_value'] for x in rows if x['group']=='control' and x['log2_value'] is not None];inf=[x['log2_value'] for x in rows if x['group']=='inflamed' and x['log2_value'] is not None]
  r=dict(endpoint=endpoint['id'],model=model,n_control=len(ctrl),n_inflamed=len(inf),evaluable=min(len(ctrl),len(inf))>=2,protein_adjusted=False)
  if r['evaluable']:
   delta=statistics.mean(inf)-statistics.mean(ctrl);median=statistics.median(ctrl);ld=[]
   for x in rows:
    if x['log2_value'] is None:continue
    a=[z['log2_value'] for z in rows if z['group']=='inflamed' and z['sample']!=x['sample'] and z['log2_value'] is not None];b=[z['log2_value'] for z in rows if z['group']=='control' and z['sample']!=x['sample'] and z['log2_value'] is not None]
    d=statistics.mean(a)-statistics.mean(b);ld.append(d);loo.append(dict(endpoint=endpoint['id'],model=model,left_out_sample=x['sample'],mean_log2_difference=d))
   r.update(mean_log2_difference=delta,geometric_mean_fold_change=2**delta,median_log2_difference=statistics.median(inf)-median,inflamed_above_control_median=sum(v>median for v in inf),inflamed_below_control_median=sum(v<median for v in inf),LOO_min=min(ld),LOO_max=max(ld))
  else:r['reason']='fewer than two evaluable animals in at least one group; no effect estimate or test'
  results.append(r)
save('RESULTS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),results=results,independent_QA_pending=True,P_values_computed=False,whole_project_complete=False))
save('ANIMAL_VALUES.json',values);save('LOO_RESULTS.json',loo)
for n,rows in [('ANIMAL_VALUES.tsv',values),('LOO_RESULTS.tsv',loo)]:
 with (W/n).open('w',encoding='utf-8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
print(json.dumps(results))
