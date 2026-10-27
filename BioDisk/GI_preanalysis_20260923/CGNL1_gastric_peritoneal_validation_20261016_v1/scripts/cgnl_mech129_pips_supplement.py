from pathlib import Path
import json,requests,zipfile,io,xml.etree.ElementTree as ET,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12518608/supplementaryFiles';a=requests.get(url,timeout=(15,60));a.raise_for_status();p=O/'literature/PIPS_supplement.zip';p.write_bytes(a.content);z=zipfile.ZipFile(io.BytesIO(a.content));assert z.testzip() is None;out={'url':url,'zip_sha256':hashlib.sha256(a.content).hexdigest(),'archive_members':z.namelist(),'documents':[]};ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
for name in z.namelist():
 if not name.endswith('.docx'):continue
 raw=z.read(name);dst=O/'literature'/Path(name).name;dst.write_bytes(raw);doc=zipfile.ZipFile(io.BytesIO(raw));assert doc.testzip() is None;root=ET.fromstring(doc.read('word/document.xml'));rows=[]
 for node in root.find('w:body',ns):
  if node.tag.endswith('}p'):rows.append({'paragraph':' '.join(t.text or '' for t in node.findall('.//w:t',ns))})
  elif node.tag.endswith('}tbl'):rows.append({'table':[[' '.join(t.text or '' for t in cell.findall('.//w:t',ns)) for cell in row.findall('w:tc',ns)] for row in node.findall('w:tr',ns)]})
 rels=doc.read('word/_rels/document.xml.rels').decode('utf-8') if 'word/_rels/document.xml.rels' in doc.namelist() else ''
 out['documents'].append({'file':str(dst),'sha256':hashlib.sha256(raw).hexdigest(),'content_in_order':rows,'relationships':rels,'media_members':[n for n in doc.namelist() if n.startswith('word/media/')],'scope':'read-only OOXML paragraph and table extraction; no pagination/layout claim'})
(O/'metadata/PIPS_SUPPLEMENT_129.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False)[:20000])
def get(acc):
 params={'accession':acc,'result':'read_run','fields':'run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format','format':'tsv'};a=requests.get('https://www.ebi.ac.uk/ena/portal/api/filereport',params=params,timeout=(15,45));a.raise_for_status();(O/'metadata'/(acc+'_submitted_files.tsv')).write_bytes(a.content);return a.text
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 for x in ex.map(get,['SRR30190686','SRR30190685']):print('SUBMITTED',x)
