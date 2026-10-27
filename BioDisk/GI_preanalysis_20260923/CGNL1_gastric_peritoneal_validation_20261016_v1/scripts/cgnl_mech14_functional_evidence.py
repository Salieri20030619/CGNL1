from pathlib import Path
import xml.etree.ElementTree as ET,json,re
O=Path.cwd()/'mechanism_investigation_v2';out={}
for fn in ['PMC10184829.xml','PMC10202830.xml','cgnl_ko_europepmc.response']:
 root=ET.fromstring((O/'literature'/fn).read_bytes());items=[]
 for sec in root.iter('sec'):
  title=sec.find('title')
  if title is not None and any(k in ''.join(title.itertext()).lower() for k in ['results','discussion']):
   for p in sec.iter('p'):
    t=' '.join(p.itertext())
    if any(x in t.lower() for x in ['rescue','rescued','barrier','stiffness','proliferation','cgnl1-ko','paracingulin-ko']):items.append(t)
 out[fn]=list(dict.fromkeys(items))
 print(fn,'\n'.join(out[fn])[:18000])
(O/'literature/FUNCTIONAL_RELEVANT_PARAGRAPHS.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
