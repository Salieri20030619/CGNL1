from pathlib import Path
import xml.etree.ElementTree as ET,json
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'metadata/HPA_CGNL1.xml').read_bytes());parents={c:p for p in root.iter() for c in p};records=[]
for e in root.iter('tissue'):
 if 'stomach' not in (e.text or '').lower():continue
 p=parents[e];path=[];cur=e
 while cur in parents:
  cur=parents[cur];path.append({'tag':cur.tag,'attributes':cur.attrib})
 rec={'tissue':e.text,'path':path,'parent_xml':ET.tostring(p,encoding='unicode')};records.append(rec)
(O/'metadata/HPA_STOMACH_PROTEIN_RECORDS.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8');print('RECORDS',len(records))
for r in records:
 if any('cancer' in str(p).lower() for p in r['path']) or 'cancer' in r['tissue'].lower():print(json.dumps(r,ensure_ascii=False)[:10000])
print('ROOT_TAG',root.tag,'ENTRY_TOP_TAGS',[(e.tag,e.attrib) for e in root.findall('.//entry/*')][:25])
