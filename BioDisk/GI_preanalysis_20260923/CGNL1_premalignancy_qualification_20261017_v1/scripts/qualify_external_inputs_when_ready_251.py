import pathlib,json,gzip,csv,subprocess
csv.field_size_limit(10000000)
p=pathlib.Path.cwd();d=json.loads((p/'metadata/EXTERNAL_DOWNLOAD_LAUNCH_250.json').read_text())
subprocess.run(['/usr/bin/python3',str(p.parent/'CGNL1_gastric_peritoneal_validation_20261016_v1/scripts/cgnl_wait_pidfd_system.py'),str(d['wrapper_pid'])],check=True)
assert json.loads(pathlib.Path(d['receipt']).read_text())['exit_code']==0
assert json.loads((p/'metadata/EXTERNAL_DOWNLOAD_COMPLETE_250.json').read_text())['all_pass']
raw=p/'raw/external_250';out=p/'processed/external_251';out.mkdir(exist_ok=True)
metadata=[];samples={};header=None;n=0
with gzip.open(raw/'GSE78523_series_matrix.txt.gz','rt',encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0].startswith('!'):
   if row[0].startswith('!Sample_'):metadata.append(row)
   continue
  if header is None:
   header=row;assert row[0]=='ID_REF' and len(set(row[1:]))==45;continue
  assert len(row)==len(header);n+=1
assert n>1000
(out/'sample_metadata_rows.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
probes=[];platform_header=None;inside=False
with gzip.open(raw/'GPL18990_family.soft.gz','rt',encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!platform_table_begin':inside=True;continue
  if row[0]=='!platform_table_end':inside=False;continue
  if not inside:continue
  if platform_header is None:platform_header=row;continue
  if any('CGNL1' in x.upper() or 'PARACINGULIN' in x.upper() or 'JACOP'==x.upper() for x in row):probes.append(row)
assert platform_header is not None
(out/'CGNL1_annotation_candidates.json').write_text(json.dumps({'header':platform_header,'rows':probes},indent=2),encoding='utf-8')
result={'qualification_only':True,'samples':len(header)-1,'expression_rows':n,'candidate_annotation_rows':len(probes),'statistics_run':False,'requires_review':'patient labels, unique mapping, batch and clinical metadata'}
(out/'QUALIFICATION_COMPLETE_251.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
