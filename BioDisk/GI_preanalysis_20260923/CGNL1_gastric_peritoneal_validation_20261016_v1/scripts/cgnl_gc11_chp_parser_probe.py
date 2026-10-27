from pathlib import Path
import struct,json,datetime,subprocess,hashlib,shutil
O=Path.cwd()
files=list((O/'raw/original31').glob('*.chp'));assert files
p=next(p for p in files if p.name.startswith('Tumor with PC'))
b=p.read_bytes();pos=0
def take(n):
 global pos
 a=b[pos:pos+n];assert len(a)==n;pos+=n;return a
def num(fmt):return struct.unpack('>'+fmt,take(struct.calcsize('>'+fmt)))[0]
def string8():return take(num('i')).decode('utf-8')
def string16():return take(num('i')*2).decode('utf-16-be')
def params():
 out=[]
 for _ in range(num('i')):
  name=string16();value=take(num('i'));typ=string16();out.append({'name':name,'type':typ,'value_hex':value.hex(),'text_utf16':value.decode('utf-16-be').rstrip('\0') if typ=='text/plain' else None})
 return out
def header():
 out={'type':string8(),'file_id':string8(),'date':string16(),'locale':string16(),'params':params()};out['parents']=[header() for _ in range(num('i'))];return out
magic=num('B');version=num('B');ng=num('i');group_pos=num('I');h=header();groups=[]
for _ in range(ng):
 pos=group_pos;next_group=num('I');first_set=num('I');nset=num('i');name=string16();sets=[];dataset_pos=first_set
 for _ in range(nset):
  pos=dataset_pos;start=num('I');nxt=num('I');name_ds=string16();pars=params();nc=num('I');cols=[]
  for _ in range(nc):cols.append({'name':string16(),'type':num('b'),'size':num('i')})
  nr=num('I');sets.append({'name':name_ds,'start':start,'next':nxt,'params':pars,'columns':cols,'rows':nr,'first_data_hex':b[start:start+100].hex()});dataset_pos=nxt
 groups.append({'name':name,'sets':sets});group_pos=next_group
out={'file':str(p),'magic':magic,'version':version,'header':h,'groups':groups};(O/'metadata/CHP_STRUCTURE.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
url='https://s3.amazonaws.com/IonTorrent/usbupdate/TAC_4.0/hg19_AmpliSeq_Transcriptome_21K_v1.transcript.csv';dest=O/'metadata/Ion_AmpliSeq_21K_annotation.csv';r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','90','-sS','-w','%{http_code}','-o',str(dest),url],capture_output=True,timeout=95);print('ANNOTATION',r.returncode,r.stdout.decode(),dest.stat().st_size if dest.exists() else 0)
(O/'metadata/ANNOTATION_DOWNLOAD.json').write_text(json.dumps({'url':url,'exit':r.returncode,'http':r.stdout.decode(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest() if dest.exists() else None},indent=2),encoding='utf-8')
rs=shutil.which('Rscript');print('RSCRIPT',rs)
if rs:
 rr=subprocess.run([rs,'--vanilla','-e','cat(requireNamespace("affxparser", quietly=TRUE))'],capture_output=True,timeout=30);print('AFFXPARSER',rr.returncode,rr.stdout.decode(),rr.stderr.decode())
