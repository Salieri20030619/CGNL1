from pathlib import Path
import requests,json,xml.etree.ElementTree as ET,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'metadata/Yonemura_Elsevier.xml').read_bytes());text=' '.join(root.itertext());(O/'literature/Yonemura_publisher_fulltext.txt').write_text(text,encoding='utf-8');r={'text_characters':len(text),'attachments':[],'relevant_paragraphs':[]}
for e in root.iter():
 tag=e.tag.split('}')[-1]
 if tag in ['link','attachment','e-component','inter-ref'] and any(k in str(e.attrib) for k in ['mmc','supp','xlsx','attachment']):r['attachments'].append({'tag':tag,'attributes':e.attrib,'text':' '.join(e.itertext())[:300]})
 if tag in ['para','simple-para']:
  t=' '.join(e.itertext())
  if any(k in t for k in ['RNA-seq','CGNL1','DRA014736','72 h','TGF-β1','Table S']):r['relevant_paragraphs'].append(t)
(O/'metadata/YONEMURA_PRIMARY_DETAILS_73.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('ATTACHMENTS',r['attachments']);print('METHODS','\n'.join(r['relevant_paragraphs'])[-12500:])
base='https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA014/DRA014736/'
def get(n):
 a=requests.get(base+n,timeout=(15,60));(O/'metadata'/('DDBJ_'+n.replace('/','_'))).write_bytes(a.content);print('DDBJ',n,a.status_code,a.content.decode('utf-8')[:5000],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(get,['DRA014736.run.xml','DRX386385/']))
