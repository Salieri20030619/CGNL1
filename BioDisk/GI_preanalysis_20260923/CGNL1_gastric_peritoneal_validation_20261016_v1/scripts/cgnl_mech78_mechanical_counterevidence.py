from pathlib import Path
import json,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'literature/PMC10202830.xml').read_bytes());out=[]
for e in root.iter('p'):
 t=' '.join(e.itertext())
 if 'CGNL1' in t and any(k in t.lower() for k in ['fragmentation','blebbistatin','stiffness was','cg nl']):out.append(t)
(O/'literature/CGNL1_FUNCTION_CONTEXT_78.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
for i,t in enumerate(out):print('PARAGRAPH',i,t[:2600])
