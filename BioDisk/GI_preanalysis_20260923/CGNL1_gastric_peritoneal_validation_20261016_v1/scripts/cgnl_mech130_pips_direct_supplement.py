from pathlib import Path
import requests,xml.etree.ElementTree as ET,json,zipfile,io,hashlib,concurrent.futures
from urllib.parse import quote
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'literature/PMC12518608.xml').read_bytes());doi=next(x.text for x in root.findall('.//article-id') if x.attrib.get('pub-id-type')=='doi');name=next(x.attrib['{http://www.w3.org/1999/xlink}href'] for x in root.iter('media') if x.attrib.get('{http://www.w3.org/1999/xlink}href','').endswith('.docx'));url='https://media.springernature.com/original/springer-static/esm/'+quote('art:'+doi,safe='')+'/MediaObjects/'+name
a=requests.get(url,timeout=(15,60));a.raise_for_status();p=O/'literature'/name;p.write_bytes(a.content);z=zipfile.ZipFile(io.BytesIO(a.content));assert z.testzip() is None;root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};out={'url':url,'file':str(p),'sha256':hashlib.sha256(a.content).hexdigest(),'content_in_order':[],'media_members':[n for n in z.namelist() if n.startswith('word/media/')],'scope':'read-only paragraph/table structure, no pagination or layout claim'}
for node in root.find('w:body',ns):
 if node.tag.endswith('}p'):out['content_in_order'].append({'paragraph':' '.join(t.text or '' for t in node.findall('.//w:t',ns))})
 elif node.tag.endswith('}tbl'):out['content_in_order'].append({'table':[[' '.join(t.text or '' for t in cell.findall('.//w:t',ns)) for cell in row.findall('w:tc',ns)] for row in node.findall('w:tr',ns)]})
(O/'metadata/PIPS_SUPPLEMENT_130.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False)[:22000])
def get(acc):
 params={'accession':acc,'result':'read_run','fields':'run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format','format':'tsv'};a=requests.get('https://www.ebi.ac.uk/ena/portal/api/filereport',params=params,timeout=(15,45));a.raise_for_status();(O/'metadata'/(acc+'_submitted_files.tsv')).write_bytes(a.content);return a.text
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 for x in ex.map(get,['SRR30190686','SRR30190685']):print('SUBMITTED',x)
