from pathlib import Path
import json,re,xml.etree.ElementTree as ET,datetime
O=Path.cwd()/'mechanism_investigation_v2';out={}
for name in ['GSE228955_samples.soft','GSE150833.soft','GSE76162.soft']:
 p=O/'raw'/name
 if not p.exists():continue
 s=p.read_text(encoding='utf-8');keep=[l for l in s.splitlines() if l.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics','!Sample_source_name','!Sample_description','!Sample_data_processing','!Series_title','!Series_summary','!Series_overall_design','!Series_supplementary'))];out[name]=keep;print(name,'\n'.join(keep))
p=O/'metadata/PMC12775079.xml';root=ET.fromstring(p.read_bytes())
secs=[]
for sec in root.iter('sec'):
 t=sec.find('title')
 if t is not None and any(x in ''.join(t.itertext()).lower() for x in ['sample','single-cell','single cell','data availability']):secs.append(' '.join(sec.itertext()))
out['GSE308231_paper_methods']=secs;print('PAPER_METHODS',secs)
out['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(O/'metadata/DESIGN_QUALIFICATION_08.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
