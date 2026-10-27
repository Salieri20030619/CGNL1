from pathlib import Path
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
source=(B/'scripts/inspect_remote_h5_523.py').read_text(encoding='utf-8')
old="""                request=urllib.request.Request(item['content_details']['download_url'],headers={'Range':f'bytes={a}-{z}'})
                with urllib.request.urlopen(request,timeout=90) as response:
                    assert response.status==206 and response.headers.get('Content-Range')==f'bytes {a}-{z}/{self.size}','Range response mismatch'
                    data=response.read(z-a+2)
                    assert len(data)==z-a+1
                    receipt=dict(start=a,end=z,etag=response.headers.get('ETag'),sha256=hashlib.sha256(data).hexdigest())"""
new="""                import subprocess
                header=cache/(str(a)+'.headers')
                response=subprocess.run(['curl','-q','--proxy','','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','90','--range',f'{a}-{z}','--dump-header',str(header),item['content_details']['download_url']],capture_output=True)
                assert response.returncode==0,response.stderr.decode('utf-8')
                headers=header.read_text(encoding='iso-8859-1')
                ranges=re.findall(r'(?im)^content-range:\\s*(.+)$',headers)
                assert ranges and ranges[-1].strip()==f'bytes {a}-{z}/{self.size}','Range response mismatch'
                data=response.stdout;assert len(data)==z-a+1
                tags=re.findall(r'(?im)^etag:\\s*(.+)$',headers)
                receipt=dict(start=a,end=z,etag=tags[-1].strip() if tags else None,sha256=hashlib.sha256(data).hexdigest())"""
assert source.count(old)==1
source=source.replace(old,new).replace('metadata_range_cache_523','metadata_range_cache_524').replace('REMOTE_H5_METADATA_523','REMOTE_H5_METADATA_524')
exec(compile(source,'metadata_probe524','exec'))

