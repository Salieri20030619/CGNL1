from pathlib import Path
import sys,json,urllib.request,zipfile,struct,hashlib,zlib,time,concurrent.futures
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024')
sys.path.insert(0,str(P))
from remote_zip_reader import Remote

def get_range(url,start,size,total):
    req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{start+size-1}'})
    with urllib.request.urlopen(req,timeout=60) as r:
        assert r.status==206
        assert r.headers['Content-Range']==f'bytes {start}-{start+size-1}/{total}'
        data=r.read()
    assert len(data)==size
    return data

def acquire(item):
    label,member=item;folder=P/label;folder.mkdir(exist_ok=True)
    dest=folder/Path(member).name
    receipt=dest.with_name(dest.name+'.receipt.json')
    if dest.exists():
        rec=json.loads(receipt.read_text());assert hashlib.sha256(dest.read_bytes()).hexdigest()==rec['sha256']
        return {'set':label,'member':member,'reused':True}
    url=f'https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/06/PXD059245/{label}.zip'
    rr=Remote(url)
    with zipfile.ZipFile(rr) as z:info=z.getinfo(member)
    header=get_range(url,info.header_offset,30,rr.size)
    assert header[:4]==b'PK\x03\x04'
    n,e=struct.unpack_from('<HH',header,26)
    start=info.header_offset+30+n+e
    spec={'url':url,'archive_bytes':rr.size,'member':member,'start':start,'compressed_bytes':info.compress_size,'crc32':info.CRC,'bytes':info.file_size,'method':info.compress_type}
    meta=dest.with_name(dest.name+'.range_spec.json')
    if meta.exists():assert json.loads(meta.read_text())==spec
    else:meta.write_text(json.dumps(spec,indent=2),encoding='utf-8')
    partial=dest.with_name(dest.name+'.deflate.partial')
    done=partial.stat().st_size if partial.exists() else 0
    assert done<=info.compress_size
    with partial.open('ab') as h:
        while done<info.compress_size:
            count=min(1024*1024,info.compress_size-done)
            for attempt in range(2):
                try:data=get_range(url,start+done,count,rr.size);break
                except Exception:
                    if attempt==1:raise
            h.write(data);h.flush();done+=len(data)
            print(json.dumps({'set':label,'member':member,'compressed_saved':done,'total':info.compress_size}),flush=True)
    compressed=partial.read_bytes()
    data=zlib.decompress(compressed,-15) if info.compress_type==8 else compressed
    assert len(data)==info.file_size and zlib.crc32(data)==info.CRC
    staging=dest.with_name(dest.name+'.verified_tmp');staging.write_bytes(data);staging.rename(dest)
    rec={**spec,'sha256':hashlib.sha256(data).hexdigest(),'compressed_sha256':hashlib.sha256(compressed).hexdigest(),'crc_pass':True}
    receipt.write_text(json.dumps(rec,indent=2),encoding='utf-8')
    return {'set':label,'member':member,'complete':True,'sha256':rec['sha256']}

if __name__=='__main__':
    results=[]
    for item in [('Set03_Phospho','txt/Phospho (STY)Sites.txt'),('Set04_Phospho','txt/evidence.txt')]:
        try:result=acquire(item)
        except Exception as e:result={'set':item[0],'member':item[1],'error':repr(e)}
        results.append(result);print(json.dumps(result),flush=True)
        (P/'CHUNKED_ACQUISITION_STATUS.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
