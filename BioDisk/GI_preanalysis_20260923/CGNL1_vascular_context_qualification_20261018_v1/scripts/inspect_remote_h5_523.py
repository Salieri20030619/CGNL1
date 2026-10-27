from pathlib import Path
import collections,datetime,hashlib,io,json,re,urllib.request,traceback
import h5py,numpy as np
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
manifest=json.loads((B/'metadata/external_sources_521/mendeley_files.json').read_text(encoding='utf-8'))
item=next(x for x in manifest if x['filename'].endswith('.h5ad'))
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/independent_processed_522')
cache=D/'metadata_range_cache_523';cache.mkdir(exist_ok=True)
prefix=D/(item['filename']+'.part')
class Remote(io.RawIOBase):
    def __init__(self):
        self.pos=0;self.size=item['size'];self.block=65536;self.receipts=[]
        self.prefix_size=prefix.stat().st_size;self.transferred=0
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
        if self.pos<0:raise ValueError('Negative offset')
        return self.pos
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n<=0:return b''
        assert n<=8*1024*1024,'Unexpected large metadata request'
        start=self.pos;end=start+n;out=[]
        if end<=self.prefix_size:
            with prefix.open('rb') as f:f.seek(start);data=f.read(n)
            assert len(data)==n;self.pos=end;return data
        while self.pos<end:
            a=self.pos//self.block*self.block;z=min(a+self.block,self.size)-1
            p=cache/(str(a)+'.bin')
            if p.exists():data=p.read_bytes()
            else:
                assert self.transferred<16*1024*1024,'Metadata transfer budget reached'
                request=urllib.request.Request(item['content_details']['download_url'],headers={'Range':f'bytes={a}-{z}'})
                with urllib.request.urlopen(request,timeout=90) as response:
                    assert response.status==206 and response.headers.get('Content-Range')==f'bytes {a}-{z}/{self.size}','Range response mismatch'
                    data=response.read(z-a+2)
                    assert len(data)==z-a+1
                    receipt=dict(start=a,end=z,etag=response.headers.get('ETag'),sha256=hashlib.sha256(data).hexdigest())
                p.write_bytes(data);self.receipts.append(receipt);self.transferred+=len(data)
            take=min(end-self.pos,len(data)-(self.pos-a));assert take>0
            out.append(data[self.pos-a:self.pos-a+take]);self.pos+=take
        return b''.join(out)
    def readinto(self,b):
        data=self.read(len(b));b[:len(data)]=data;return len(data)
def decode(v):return v.decode('utf-8') if isinstance(v,bytes) else v.item() if isinstance(v,np.generic) else v
remote=Remote();result=dict(scope='Metadata-only inspection; whole-file SHA256 not yet validated',source_manifest=item,scientific_validation=False)
try:
    with h5py.File(remote,'r') as f:
        result['root_keys']=list(f);result['obs_fields']=list(f['obs'])
        result['obs']={}
        for key,obj in f['obs'].items():
            if not re.search('patient|sample|cell.?type|hgp|growth|treat|annot|orig|batch|donor|cluster',key,re.I):continue
            if isinstance(obj,h5py.Group) and 'categories' in obj and 'codes' in obj:
                cats=[decode(v) for v in obj['categories'][()]]
                counts=collections.Counter(obj['codes'][()].tolist())
                result['obs'][key]={str(cats[i]) if i>=0 else '<missing>':int(n) for i,n in counts.items()}
            elif isinstance(obj,h5py.Dataset):
                values=[decode(v) for v in obj[()]];counts=collections.Counter(values)
                result['obs'][key]=dict(unique=len(counts),counts=dict(counts) if len(counts)<=100 else None)
        result['var_fields']=list(f['var'])
        result['layers']=list(f['layers']) if 'layers' in f else []
        result['X_attributes']={k:decode(v) if np.ndim(v)==0 else [decode(i) for i in v] for k,v in f['X'].attrs.items()}
        result['status']='metadata_read'
except Exception:result.update(status='failed',error=traceback.format_exc())
result.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),new_downloaded_bytes=remote.transferred,range_receipts=remote.receipts,prefix_bytes_available=remote.prefix_size)
(B/'metadata/REMOTE_H5_METADATA_523.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ('source_manifest','range_receipts')}),flush=True)

