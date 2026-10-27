from pathlib import Path
import gzip,collections,json,hashlib
O=Path.cwd()/'mechanism_investigation_v2';out={}
for p in sorted((O/'raw').glob('GSE223165*.gz')):
 try:
  b=gzip.decompress(p.read_bytes());t=b.decode('utf-8');ls=t.splitlines();counts=collections.Counter(len(l.split('\t')) for l in ls);expected=len(ls[0].split('\t'));bad=[{'line':i+1,'fields':len(l.split('\t')),'prefix':l[:260],'suffix':l[-140:]} for i,l in enumerate(ls) if len(l.split('\t'))!=expected]
  v={'compressed_bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'uncompressed_bytes':len(b),'header_fields':expected,'line_widths':dict(counts),'bad_rows':bad,'gzip_crc':True,'null_bytes':b.count(b'\x00')};out[p.name]=v;print(p.name,'WIDTHS',dict(counts),'BAD',bad[:4]);print('TARGET_ROWS',[l[:650] for l in ls if l.startswith('ENSG00000128849\t')])
 except Exception as e:out[p.name]={'error':str(e)};print(p.name,str(e))
(O/'metadata/EXPOSURE_INPUT_INTEGRITY_25.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
