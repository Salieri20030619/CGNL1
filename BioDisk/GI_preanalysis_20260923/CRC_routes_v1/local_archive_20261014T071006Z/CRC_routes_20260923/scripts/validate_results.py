import pathlib,json,csv,gzip,statistics,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8');O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');checks=[]
def ck(name,v):assert v,name;checks.append({'check':name,'passed':True})
def read(p):
 f=gzip.open if str(p).endswith('.gz') else open
 with f(p,'rt',encoding='utf-8-sig',newline='') as q:return list(csv.DictReader(q,delimiter='\t'))
for r in ['route1','route2']:ck(r+' exit marker',json.loads((O/'results'/r/'COMPLETE.json').read_text())['exit_code']==0)
ck('integration exit',json.loads((O/'results/INTEGRATION_COMPLETE.json').read_text())['exit_code']==0)
for route,keys in [('route1',['celltype','program']),('route2',['contrast','program'])]:
 p=read(O/'results'/route/'program_patient_effects.tsv');s=read(O/'results'/route/'program_summary.tsv')
 for row in s:
  vals=[float(x['delta']) for x in p if all(x[k]==row[k] for k in keys) and x['delta'] not in ['', 'nan']]
  if not vals:ck(route+' unavailable proxy preserved',row['median_delta'] in ['','nan']);continue
  ck(route+' patient median '+str([row[k] for k in keys]),abs(statistics.median(vals)-float(row['median_delta']))<1e-10)
meta=read(O/'results/route2/sample_QC.tsv');ck('8 sections and5 patients',len(meta)==8 and len(set(x['patient'] for x in meta))==5);ck('60801 matched annotated spots',sum(int(x['annotated_spots']) for x in meta)==60801)
el=read(O/'results/route2/contrast_eligibility.tsv')
for contrast,pts,secs in [('tumor_front_vs_other_tumor',5,7),('tumor_front_vs_deep_tumor',4,6)]:
 z=[r for r in el if r['contrast']==contrast and r['eligible']=='True'];ck(contrast+' count',len(z)==secs and len({r['patient'] for r in z})==pts)
a=read(O/'results/route1/cross_subtype_sensitivity_screen.tsv');passes=[r for r in a if r['passes_descriptive_sensitivity']=='True'];ck('126 descriptive sensitivity candidates',len(passes)==126);ck('strictQC effect threshold',all(abs(float(r['Tumor_strict_QC_delta']))>=1 for r in passes))
records=[]
for name,url in json.loads((O/'spatial_download_plan.json').read_text()):
 choices=[]
 for suffix in ['.complete.json','.retry.complete.json','.openssl.complete.json']:
  f=O/'raw'/(name+suffix)
  if f.exists():
   j=json.loads(f.read_text());choices.extend(j if isinstance(j,list) else [j])
 good=[r for r in choices if r.get('exit_code')==0 and r.get('transfer_exit_code')==0];ck('download '+name,bool(good));r=good[-1];r['file']=name+'.tar.gz';records.append(r)
j=json.loads((O/'raw/COUNTS_DOWNLOAD.json').read_text());j['file']='nonimmune_counts.txt.gz';records.append(j);ck('9 input files complete',len(records)==9)
(O/'FINAL_INPUT_MANIFEST.json').write_text(json.dumps(records,indent=2));(O/'QA.json').write_text(json.dumps({'checks':checks,'passed':len(checks),'visual_QA':'Three figures and original clinical PDF table page13 visually inspected; no plotting inference tests.'},ensure_ascii=False,indent=2),encoding='utf-8')
print('PASSED',len(checks));print('DATA_GB',sum(x['bytes'] for x in records)/1e9);print('ROUTE2_SCREEN',json.loads((O/'results/INTEGRATION_COMPLETE.json').read_text()));print('TU_COVERAGE',[r for r in read(O/'results/route1/paired_coverage.tsv') if r['celltype'] in ['Tumor_author','Tumor_strict_QC','Tu01_AREG','Tu04_RGMB']])
