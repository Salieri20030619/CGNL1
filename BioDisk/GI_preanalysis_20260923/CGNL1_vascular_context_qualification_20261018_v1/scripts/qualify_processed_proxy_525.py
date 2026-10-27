from pathlib import Path
import json,os,datetime
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/independent_processed_522')
part=D/'CRCLM_multiome_GEX_decontaminated.h5ad.part'
assert part.exists() and not part.with_suffix('').exists()
with (B/'metadata/PROXY525_LAUNCH.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),started=datetime.datetime.now(datetime.timezone.utc).isoformat(),resume_bytes=part.stat().st_size,route='sgs10 reverse SSH tunnel to authorized Windows proxy',partial_preserved=True),f,indent=2)
source=(B/'scripts/qualify_processed_522.py').read_text(encoding='utf-8')
assert source.count("'--proxy','','--noproxy','*'")==1
source=source.replace("'--proxy','','--noproxy','*'","'--proxy','http://127.0.0.1:7897','--noproxy',''")
source=source.replace('D.mkdir(exist_ok=False)','D.mkdir(exist_ok=True)')
source=source.replace('PROCESSED522_STATE','PROCESSED525_STATE').replace("D/'download.log'","D/'download_proxy525.log'")
exec(compile(source,'processed525_proxy_resume','exec'))

