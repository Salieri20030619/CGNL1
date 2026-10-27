from pathlib import Path
import requests,json,csv,gzip,io,hashlib,threading,concurrent.futures,datetime,time
p=Path.cwd();out=p/'raw/background_check_312';out.mkdir(exist_ok=True);meta={r['accession']:r for r in json.loads((p/'metadata/GSE140973_SAMPLE_METADATA_295.json').read_text())}
with (p/'results/paired_pdx_310/target_background_QC.tsv').open(encoding='utf-8') as f:tasks=[r for r in csv.DictReader(f,delimiter='\t') if r['above_negative99'].lower()!='true']
assert len(tasks)==14;tls=threading.local()
def fetch(task):
 mid=task['accession'];url=meta[mid]['!Sample_supplementary_file'][0].replace('ftp://','https://',1);assert url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/samples/')
 if not hasattr(tls,'s'):tls.s=requests.Session();tls.s.trust_env=False
 error=None
 for attempt in range(1,3):
  try:
   r=tls.s.get(url,timeout=(15,90));r.raise_for_status();raw=gzip.decompress(r.content)
   try:text=raw.decode('utf-8');encoding='utf-8'
   except UnicodeDecodeError:text=raw.decode('cp1252');encoding='cp1252';assert text.encode(encoding)==raw
   f=out/url.rsplit('/',1)[1];assert not f.exists();f.write_bytes(r.content);header=None;found=[];params=None;protocol=None
   for row in csv.reader(io.StringIO(text),delimiter='\t'):
    if not row:continue
    if row[0]=='FEPARAMS':params=row;continue
    if params is not None and protocol is None and row[0]=='DATA':protocol=dict(zip(params,row))
    if row[0]=='FEATURES':header=row;continue
    if header and row[0]=='DATA':
     d=dict(zip(header,row))
     if d['FeatureNum']=='22249':found.append(d)
   assert len(found)==1 and found[0]['ProbeName']=='A_23_P163306';d=found[0];assert float(d['gProcessedSignal'])==float(task['CGNL1_raw_signal']);assert protocol['Scan_NumChannels']=='1' and 'gProcessedSignal' in header and 'rProcessedSignal' not in header
   return dict(task,success=True,url=url,path=str(f),sha256=hashlib.sha256(r.content).hexdigest(),gzip_CRC_pass=True,encoding=encoding,protocol=protocol['Protocol_Name'],**{k:d[k] for k in ['gIsWellAboveBG','gIsPosAndSignif','gIsFeatNonUnifOL','gIsSaturated','gBGMeanSignal']})
  except Exception as exc:
   error=repr(exc)
   if isinstance(exc,requests.HTTPError) and exc.response.status_code in [401,403]:break
   if attempt<2:time.sleep(15)
 return dict(task,success=False,error=error)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(fetch,tasks))
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selection':'All14arrays flagged by prespecified negative-control99percentile QC, independent of biological effect/P value','records':records,'all_pass':all(r['success'] for r in records),'vendor_well_above_background':sum(r.get('gIsWellAboveBG')=='1' for r in records),'nonuniform_flags':sum(r.get('gIsFeatNonUnifOL')=='1' for r in records),'saturated_flags':sum(r.get('gIsSaturated')=='1' for r in records),'whole_project_complete':False};(p/'results/paired_pdx_310/RAW_BACKGROUND_QA_312.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in x.items() if k!='records'}),flush=True);assert x['all_pass']
