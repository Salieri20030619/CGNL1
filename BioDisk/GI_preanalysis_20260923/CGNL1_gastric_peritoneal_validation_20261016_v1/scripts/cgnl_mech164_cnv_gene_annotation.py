from pathlib import Path
import requests,json,hashlib,gzip,datetime,concurrent.futures,re,threading,time
O=Path.cwd()/'mechanism_investigation_v2';R=O/'reference/gencode_v47';name='gencode.v47.primary_assembly.basic.annotation.gtf.gz';md={x.split()[1].lstrip('*'):x.split()[0] for x in (R/'MD5SUMS').read_text().splitlines() if len(x.split())==2};expected=md[name];url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/'+name
with requests.get(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},stream=True,timeout=(15,45)) as a:assert a.status_code==206;size=int(a.headers['Content-Range'].split('/')[1])
folder=R/'cnv_annotation_ranges';folder.mkdir(exist_ok=True);manifest={};lock=threading.Lock()
def part(iv):
 start,end=iv;p=folder/f'{start}_{end}.bin'
 for attempt in range(1,4):
  tmp=p.with_name(p.name+f'.attempt{attempt}.partial')
  try:
   begin=time.monotonic();n=0
   with requests.get(url,headers={'Range':f'bytes={start}-{end-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,90)) as a:
    assert a.status_code==206 and a.headers['Content-Range']==f'bytes {start}-{end-1}/{size}'
    with tmp.open('wb') as f:
     for b in a.iter_content(65536):
      f.write(b);n+=len(b)
      if n>end-start or time.monotonic()-begin>1800:raise TimeoutError('bounded chunk')
   assert n==end-start;tmp.rename(p)
   with lock:manifest[str(start)]={'bytes':n,'success':True};(folder/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
   return p
  except Exception as e:err=repr(e)
 raise RuntimeError(err)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:parts=list(ex.map(part,[(x,min(x+8388608,size)) for x in range(0,size,8388608)]))
p=R/name;tmp=p.with_name(p.name+'.assembled.partial');h=hashlib.md5();n=0
with tmp.open('wb') as f:
 for q in parts:
  with q.open('rb') as src:
   while b:=src.read(4194304):f.write(b);h.update(b);n+=len(b)
assert n==size and h.hexdigest()==expected;assert not p.exists();tmp.rename(p);rows=[];seen=set();chroms={'chr'+str(i) for i in range(1,23)}
with gzip.open(p,'rt',encoding='utf-8') as f:
 for line in f:
  if line.startswith('#'):continue
  fields=line.rstrip('\n').split('\t')
  if fields[2]!='gene' or fields[0] not in chroms:continue
  attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]));gid=attrs['gene_id'].split('.')[0];assert gid not in seen;seen.add(gid);rows.append((gid,fields[0],int(fields[3]),int(fields[4]),attrs.get('gene_name',''),attrs.get('gene_type','')))
out=R/'autosomal_gene_order_v47.tsv'
with out.open('w',encoding='utf-8',newline='') as f:
 import csv
 w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','chromosome','start','end','gene_symbol','gene_type']);w.writerows(rows)
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':url,'bytes':size,'official_md5':expected,'MD5_pass':True,'autosomal_genes':len(rows),'gene_order':str(out),'gene_order_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'assembly':'GRCh38,GENCODE47 primaryassemblybasic;1-basedinclusiveGTFcoordinates;onlychr1-22'};(O/'metadata/CNV_ANNOTATION_READY_164.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt),flush=True)
