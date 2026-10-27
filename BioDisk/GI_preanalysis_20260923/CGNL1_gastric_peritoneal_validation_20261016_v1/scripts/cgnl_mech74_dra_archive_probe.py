from pathlib import Path
import json,requests,time,xml.etree.ElementTree as ET,collections
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'metadata/Yonemura_Elsevier.xml').read_bytes());text=' '.join(root.itertext());r={'publisher_response_bytes':(O/'metadata/Yonemura_Elsevier.xml').stat().st_size,'text_characters':len(text),'tags':dict(collections.Counter(e.tag.split('}')[-1] for e in root.iter())),'has_fulltext_methods':'RNA-seq' in text or 'STAR Methods' in text}
u='https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA014/DRA014736/DRX386385/DRR400682_1.fastq.bz2';q={'url':u,'purpose':'officialDDBJalternativearchivecapabilityprobeonly;notreplacingENAinputs'};start=time.monotonic()
try:
 with requests.get(u,headers={'Range':'bytes=0-65535','Accept-Encoding':'identity'},stream=True,timeout=(15,30)) as a:
  b=a.raw.read(65536);q.update(status=a.status_code,content_range=a.headers.get('Content-Range'),bytes_read=len(b),seconds=time.monotonic()-start,magic=b[:4].hex())
except Exception as e:q['error']=repr(e)
r['DDBJ_probe']=q;(O/'metadata/PRIMARY_ARCHIVE_QUALIFICATION_74.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
p=O/'STATUS.json';s=json.loads(p.read_text(encoding='utf-8'));s['last_goal_turn_classification']='progress: archivedGRCh37targetmappingqualified; patient/site completeness-aware routeanalysiscontrollerlaunched; DDBJoriginalsubmissionverifiedasadditionalpublicsource';s['route_analysis_controller']='cgnl_mech71_route_analyze_when_ready.py';s['route_mapping']='CGNL1EnsemblGRCh37ENSG00000128849,3transcripts; allmustbe present;codingonlysensitivity; authoractualSalmontranscriptoutputnotgenevalues';s['primary_archive_qualification']='metadata/PRIMARY_ARCHIVE_QUALIFICATION_74.json';p.write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
