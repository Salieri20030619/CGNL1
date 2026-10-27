from pathlib import Path

# Preserve the previous launcher and its receipts; reuse its checked resume logic.
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
source=B/'scripts/download_proxy_493.py'
code=source.read_text(encoding='utf-8')
old="['curl','--proxy','http://127.0.0.1:7897','--noproxy',''"
new="['curl','-q','--proxy','','--noproxy','*'"
assert code.count(old)==1
code=code.replace(old,new)
code=code.replace('PROXY493_','DIRECT494_')
code=code.replace("'.resume493.log'","'.direct494.log'")
code=code.replace("proxy='http://127.0.0.1:7897'","proxy=None,route='server direct HTTPS'")
assert 'http://127.0.0.1:7897' not in code
exec(compile(code,str(source)+' [direct494]', 'exec'))

