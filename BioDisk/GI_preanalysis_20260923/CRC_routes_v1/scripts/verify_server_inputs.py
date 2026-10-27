import pathlib,json,hashlib
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');r=json.loads((O/'metadata/FINAL_INPUT_MANIFEST.json').read_text());out=[]
for x in r:
 p=O/'raw'/x['file'];h=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();assert h==x['sha256'];assert p.stat().st_size==x['bytes'];out.append({'file':x['file'],'bytes':p.stat().st_size,'sha256':h,'matches_local':True})
(O/'results/SERVER_INPUT_INTEGRITY.json').write_text(json.dumps({'exit_code':0,'files':out},indent=2));print('All',len(out),'input SHA256 values match local download records')
