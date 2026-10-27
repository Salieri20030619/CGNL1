from pathlib import Path
import json,hashlib,stat,datetime
B=Path(r'D:\CGNL1');H=B/'handoff';U=Path(r'C:\Users\31937').resolve()
manifest=json.loads((H/'LOCAL_ARCHIVE_MANIFEST.json').read_text(encoding='utf-8'))
verified=json.loads((H/'ARCHIVE_VERIFICATION.json').read_text(encoding='utf-8'));assert verified['all_pass'] and verified['files']==len(manifest)
previous=json.loads((H/'MIGRATION_STATUS.json').read_text(encoding='utf-8'))
(H/'MIGRATION_FIRST_ATTEMPT.json').write_text(json.dumps(previous,indent=2),encoding='utf-8')
roots=sorted({U/Path(x['archive_relative']).parts[0] for x in manifest})
for root in roots:assert root.resolve().parent==U and (root.name.startswith('CGNL1') or root.name=='CRC_routes_20260923')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
remaining=[]
for x in manifest:
 p=Path(x['local']).resolve()
 assert any(r.resolve() in p.parents for r in roots)
 if p.exists():assert p.is_file() and digest(p)==x['sha256'];remaining.append(p)
readonly=[]
for p in remaining:
 try:p.unlink()
 except PermissionError:
  p.chmod(stat.S_IREAD|stat.S_IWRITE);p.unlink();readonly.append(str(p))
for root in roots:
 if not root.exists():continue
 for p in sorted([p for p in root.rglob('*') if p.is_dir()],key=lambda p:len(p.parts),reverse=True):p.rmdir()
 root.rmdir()
raw=json.loads((H/'RAW_TRANSFER_MANIFEST.json').read_text(encoding='utf-8'))
assert all(not Path(x['local']).exists() for x in manifest+raw)
receipt={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_source_files_verified_before_delete':True,'directories_removed':[str(r) for r in roots],'files_removed':len(manifest)+len(raw),'bytes_removed':sum(x['bytes'] for x in manifest+raw),'remaining_source_files':0,'readonly_attributes_cleared_after_verification':readonly,'server_archive':verified['archive']}
(H/'LOCAL_CLEANUP_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
final={'utc':receipt['completed_utc'],'stage':'complete','archive':verified['archive'],'analysis':previous['analysis'],'files_archived':receipt['files_removed'],'bytes_archived':receipt['bytes_removed'],'all_server_hashes_verified':True,'all_local_sources_removed':True,'local_final_reports':len(json.loads((H/'FINAL_REPORTS_RETAINED.json').read_text(encoding='utf-8')))}
(H/'MIGRATION_STATUS.json').write_text(json.dumps(final,indent=2),encoding='utf-8')
(H/'CLEANUP_RETRY_EXIT.json').write_text(json.dumps({'exit_code':0,'completed_utc':receipt['completed_utc']},indent=2),encoding='utf-8')
print(json.dumps(final,indent=2),flush=True)