from pathlib import Path
import struct,json,datetime,hashlib
import numpy as np,pandas as pd
O=Path.cwd();receipts=json.loads((O/'metadata/ORIGINAL31_DOWNLOAD_RECEIPTS.json').read_text());assert len(receipts)==54 and all(r['hash_matches'] for r in receipts)
class Reader:
 def __init__(self,b):self.b=b;self.p=0
 def get(self,n):a=self.b[self.p:self.p+n];assert len(a)==n;self.p+=n;return a
 def num(self,f):return struct.unpack('>'+f,self.get(struct.calcsize('>'+f)))[0]
 def s8(self):return self.get(self.num('i')).decode('utf-8')
 def sw(self):return self.get(self.num('i')*2).decode('utf-16-be')
 def params(self):
  out={}
  for _ in range(self.num('i')):
   k=self.sw();v=self.get(self.num('i'));t=self.sw();out[k]={'type':t,'value':v.decode('ascii') if t=='text/ascii' else v.decode('utf-16-be').rstrip('\0') if t=='text/plain' else v.hex()}
  return out
 def header(self):
  h={'type':self.s8(),'id':self.s8(),'date':self.sw(),'locale':self.sw(),'params':self.params()};h['parents']=[self.header() for _ in range(self.num('i'))];return h
 def parse(self):
  assert self.num('B')==59 and self.num('B')==1;ng=self.num('i');gp=self.num('I');h=self.header();assert ng==1;self.p=gp;nextg=self.num('I');dp=self.num('I');nds=self.num('i');gn=self.sw();assert nds==1;self.p=dp;start=self.num('I');nextd=self.num('I');dn=self.sw();pars=self.params();cols=[]
  for _ in range(self.num('I')):cols.append((self.sw(),self.num('b'),self.num('i')))
  n=self.num('I');assert cols[0][1]==7 and cols[1][1]==6;rowwidth=sum(c[2] for c in cols);assert rowwidth==25;names=[];vals=[]
  for i in range(n):
   self.p=start+i*rowwidth;length=self.num('i');block=self.get(cols[0][2]-4);assert length<=len(block) and all(v==0 for v in block[length:]);names.append(block[:length].decode('ascii'));vals.append(self.num('f'))
  assert len(set(names))==n;arr=np.frombuffer(self.b,dtype=[('length','>i4'),('symbol','S17'),('value','>f4')],count=n,offset=start);assert np.array_equal(arr['value'],vals);assert [s.decode('ascii') for s in arr['symbol']]==names
  assert start+n*rowwidth==nextd and nextd==len(self.b);return h,names,np.array(vals,dtype=float)
mat=[];meta=[];headers=[];gene_ref=None
for r in receipts:
 p=O/'raw/original31'/r['filename'];b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==r['expected_sha256'];h,g,y=Reader(b).parse();assert np.isfinite(y).all();gene_ref=g if gene_ref is None else gene_ref;assert g==gene_ref;mat.append(y);group='with_PC' if p.name.startswith('Tumor with PC ') else 'without_PC' if p.name.startswith('Tumor without PC ') else 'healthy_unpaired';v=lambda key:h['params'].get(key,{}).get('value');meta.append({'filename':p.name,'group_from_filename':group,'folder_name':r['folder_name'],'file_uuid':h['id'],'algorithm':v('affymetrix-algorithm-name'),'declared_scale':v('affymetrix-algorithm-param-quantification-scale'),'declared_type':v('affymetrix-algorithm-param-quantification-type'),'program_version':v('program-version'),'min':float(y.min()),'median':float(np.median(y)),'max':float(y.max()),'sum':float(y.sum()),'zeros':int((y==0).sum()),'expression_sha256':hashlib.sha256(y.tobytes()).hexdigest()});headers.append({'filename':p.name,'header':h})
m=pd.DataFrame(meta);X=pd.DataFrame(np.column_stack(mat),index=gene_ref,columns=m.filename);X.index.name='symbol';X.to_csv(O/'processed/original31_vendor_quantification.tsv.gz',sep='\t');m.to_csv(O/'metadata/original31_sample_manifest.tsv',sep='\t',index=False);(O/'metadata/ALL_CHP_HEADERS.json').write_text(json.dumps(headers,indent=2),encoding='utf-8')
qualification={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':len(m),'features':len(gene_ref),'CGNL1_present':'CGNL1' in gene_ref,'parallel_binary_decoders_agree':True,'duplicate_UUIDs':int(m.file_uuid.duplicated().sum()),'duplicate_expression_vectors':int(m.expression_sha256.duplicated().sum()),'group_sizes':m.group_from_filename.value_counts().to_dict(),'scale_range':m[['min','median','max','sum','zeros']].describe().to_dict(),'group_assignment':'filename groups match published12/19 counts;folder name conflicts;headers contain processing metadata but no clinicalgroup;clinicalindividualmapunverified','clinical_label_qualification':'provisional filename-based reproduction;retain ambiguity, do not resolve by gene outcome','target_effects_not_yet_calculated':True};(O/'metadata/CHP_DECODE_QUALIFICATION.json').write_text(json.dumps(qualification,indent=2),encoding='utf-8');print(m.groupby('group_from_filename')[['min','median','max','sum','zeros']].agg(['min','max']).to_string());print(json.dumps(qualification))
