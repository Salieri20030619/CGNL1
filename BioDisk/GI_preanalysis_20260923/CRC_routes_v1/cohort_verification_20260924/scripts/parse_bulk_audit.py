from pathlib import Path
import gzip,csv,json,re,urllib.request,numpy as np,pandas as pd,hashlib
D=Path(__file__).resolve().parent.parent/'bulk'
# Retain only the platform annotation table from the large family archive; sample records are not needed.
src=D/'GPL17586_family.soft.gz';out=D/'GPL17586_platform_table.tsv.gz'
if not out.exists():
 with gzip.open(src,'rt',encoding='utf-8') as r,gzip.open(out,'wt',encoding='utf-8') as w:
  inside=False
  for line in r:
   if line.startswith('!platform_table_begin'):inside=True;continue
   if line.startswith('!platform_table_end'):break
   if inside:w.write(line)
print('Platform table extracted',out.stat().st_size,flush=True)
p=D/'GPL6947.annot.gz';url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz'
if not p.exists():
 op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
 with op.open(url,timeout=90) as r:
  size=r.headers.get('Content-Length')
  if size:assert int(size)<100*1024*1024
  data=r.read(100*1024*1024+1);assert len(data)<100*1024*1024;p.write_bytes(data)
with gzip.open(out,'rt',encoding='utf-8') as f:
 annot=csv.DictReader(f,delimiter='\t');a=[r for r in annot if ' // CGNL1 //' in r.get('gene_assignment','')]
with gzip.open(p,'rt',encoding='utf-8') as f:
 for line in f:
  if line.startswith('ID\t'):header=line.rstrip().split('\t');break
 annot=csv.DictReader(f,fieldnames=header,delimiter='\t');b=[r for r in annot if r.get('Gene symbol')=='CGNL1']
mapping={'GSE159216':[r['ID'] for r in a],'GSE73255':[r['ID'] for r in b]};(D/'CGNL1_PLATFORM_MAPPING.json').write_text(json.dumps({'GPL17586':a,'GPL6947':b},indent=2))
summary=[]
for acc in ['GSE159216','GSE73255']:
 meta={};char=[];expr={};probes=mapping[acc];assert probes
 with gzip.open(D/(acc+'_series_matrix.txt.gz'),'rt',encoding='utf-8') as f:
  for line in f:
   if line.startswith('!Sample_'):
    z=next(csv.reader([line],delimiter='\t'));key=z[0]
    if key=='!Sample_characteristics_ch1':char.append(z[1:])
    else:meta[key]=z[1:]
   elif line.startswith('!series_matrix_table_begin'):
    ids=next(csv.reader([next(f)],delimiter='\t'))[1:]
   elif line.startswith('"') and any(line.startswith('"'+probe+'"\t') for probe in probes):
    z=next(csv.reader([line],delimiter='\t'));expr[z[0]]=np.array([float(v) if v not in ['null','NA',''] else np.nan for v in z[1:]])
 x=pd.DataFrame({'GSM':meta['!Sample_geo_accession'],'title':meta['!Sample_title']});assert list(x.GSM)==ids
 for values in char:
  if not values:continue
  keys=[v.split(': ',1)[0] for v in values if ': ' in v];key=keys[0] if keys else 'unknown';x[key]=[v.split(': ',1)[1] if ': ' in v else v for v in values]
 for probe,values in expr.items():assert len(values)==len(x);x['CGNL1__'+probe]=values
 x.to_csv(D/(acc+'_sample_metadata_CGNL1.tsv'),sep='\t',index=False)
 donorfield='patient number' if acc=='GSE159216' else 'tissue id';n=x[donorfield].nunique();missing={'NA','','--','not available','Not available','nan'}
 clinic={}
 for col in x.columns:
  if col in ['GSM','title'] or col.startswith('CGNL1__'):continue
  valid=~x[col].astype(str).isin(missing);clinic[col]={'samples_nonmissing':int(valid.sum()),'donors_nonmissing':int(x.loc[valid,donorfield].nunique()),'example_values':x.loc[valid,col].astype(str).unique()[:5].tolist()}
 rr={'dataset':acc,'expression_samples':len(x),'donor_field':donorfield,'unique_donors':n,'platforms':sorted(set(meta.get('!Sample_platform_id',[]))),'CGNL1_annotated_probes':probes,'CGNL1_measured_probes':list(expr),'CGNL1_finite_samples':{k:int(np.isfinite(v).sum()) for k,v in expr.items()},'CGNL1_range':{k:[float(np.nanmin(v)),float(np.nanmax(v))] for k,v in expr.items()},'clinical_fields':clinic,'processing':sorted(set(meta.get('!Sample_data_processing',[])))}
 summary.append(rr);print(json.dumps(rr,ensure_ascii=False),flush=True)
(D/'BULK_VERIFICATION.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'ANNOTATION_EXTRACTION.json').write_text(json.dumps({'source_file':src.name,'source_sha256':next(x['sha256'] for x in json.loads((D/'MANIFEST.json').read_text()) if x['file']==src.name),'retained_file':out.name,'retained_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'retained_bytes':out.stat().st_size,'GPL6947_url':url,'GPL6947_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'purpose':'platform annotation extracted; unrelated GEO family sample records not needed'},indent=2))
