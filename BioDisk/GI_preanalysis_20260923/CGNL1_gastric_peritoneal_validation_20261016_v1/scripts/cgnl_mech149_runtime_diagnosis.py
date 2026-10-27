from pathlib import Path
import json,datetime,time
O=Path.cwd()/'mechanism_investigation_v2';out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'ranges':[],'markers':{}}
for p in [O/'raw/GSE308231_ranges']+list((O/'raw/verified_ranges').iterdir()):
 if not p.is_dir():continue
 bins=list(p.glob('*.bin'));parts=list(p.glob('*.partial'));m=p/'range_manifest.json';q=json.loads(m.read_text()) if m.exists() else {};out['ranges'].append({'name':p.name,'verified_chunks':len(bins),'verified_bytes':sum(x.stat().st_size for x in bins),'partial_files':len(parts),'partial_bytes':sum(x.stat().st_size for x in parts),'newest_partial_age_seconds':min([time.time()-x.stat().st_mtime for x in parts],default=None),'failed_chunks':[x for x in q.values() if not x.get('success')],'newest_verified_age_seconds':min([time.time()-x.stat().st_mtime for x in bins],default=None)})
for p in [O/'metadata/GSE308231_DOWNLOAD.json',O/'results/GSE308231_SOURCE_COMPLETE.json',O/'results/GSE308231_SOURCE_QA_119.json',O/'metadata/DRA_DOWNLOAD_RECEIPTS.json',O/'metadata/DRA_QUANTIFICATION_STATUS.json']:
 if p.exists():out['markers'][p.name]=json.loads(p.read_text())
(O/'metadata/RUNTIME_DIAGNOSIS_149.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out)[:14000])
