from pathlib import Path
import json,hashlib,gzip,datetime,collections
O=Path.cwd()/'mechanism_investigation_v2';complete=json.loads((O/'metadata/ROUTE_HTTP_COMPLETE.json').read_text());rs=json.loads((O/'metadata/ROUTE_HTTP_RECEIPTS.json').read_text());assert len(rs)==66 and len({r['accession'] for r in rs})==66;verified=[];failed=[];coverage={}
for r in rs:
 patient=r['patient'];site='primary' if 'primary' in r['title'] else ('peritoneal' if 'perit' in r['title'] else ('liver' if 'liver' in r['title'] else ('lung' if 'lung' in r['title'] else ('ovary' if 'ovary' in r['title'] else 'unknown'))));key=(patient,site);coverage.setdefault(key,{'patient':patient,'site':site,'expected':0,'available':0});coverage[key]['expected']+=1
 if not r['success']:failed.append({'sample':r['accession'],'patient':patient,'site':site,'error':r.get('error')});continue
 p=O/'raw/GSE237876'/r['url'].rsplit('/',1)[1];b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==r['sha256'];raw=gzip.decompress(b);assert raw.startswith(b'Name\tLength\tEffectiveLength\tTPM\tNumReads');verified.append(r['accession']);coverage[key]['available']+=1
assert len(verified)==complete['passed'];pairs=[]
for patient in sorted(set(r['patient'] for r in rs)):
 base=coverage.get((patient,'primary'))
 if not base:continue
 for (pid,site),c in coverage.items():
  if pid==patient and site!='primary':pairs.append({'patient':pid,'site':site,'complete_pair':base['available']==base['expected'] and c['available']==c['expected'],'any_pair':base['available']>0 and c['available']>0,'primary_available':base['available'],'primary_expected':base['expected'],'met_available':c['available'],'met_expected':c['expected']})
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'download_completion':complete,'exit_success_is_not_all_files_success':True,'successful_file_hash_and_gzip_crc_verified':len(verified),'expected_files':66,'coverage':list(coverage.values()),'pairs':pairs,'failed_files':failed};(O/'metadata/ROUTE_DOWNLOAD_AUDIT_96.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print('DOWNLOAD',complete);print('PAIR_QUALIFICATION',json.dumps(pairs));print('FAILED',len(failed))
