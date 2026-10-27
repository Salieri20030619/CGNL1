from pathlib import Path
import requests,json,hashlib,zipfile,datetime
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1');D=R/'perturbation_inputs';name='GSE237653_2023_Teichmann_NCoR1_RNAseq_NormCounts.xlsx';p=D/name;partial=D/(name+'.part');s=requests.Session();s.trust_env=False;url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237653/suppl/'+name;start=partial.stat().st_size;record={'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resume_offset':start}
try:
 q=s.get(url,headers={'Range':f'bytes={start}-'},timeout=(15,90),stream=True);q.raise_for_status();record.update(status=q.status_code,content_range=q.headers.get('Content-Range'))
 if q.status_code==206:
  assert q.headers['Content-Range'].startswith(f'bytes {start}-');dest=partial;mode='ab'
 else:dest=D/(name+'.attempt2.part');mode='wb'
 with dest.open(mode) as f:
  for chunk in q.iter_content(1024*1024):f.write(chunk)
 with zipfile.ZipFile(dest) as z:assert z.testzip() is None
 dest.rename(p);record.update(integrity=True,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
except Exception as e:record['error']=repr(e)
(D/'NCOR1_RETRY_RECORD.json').write_text(json.dumps(record,indent=2));print(json.dumps(record));raise SystemExit(0 if record.get('integrity') else 1)
