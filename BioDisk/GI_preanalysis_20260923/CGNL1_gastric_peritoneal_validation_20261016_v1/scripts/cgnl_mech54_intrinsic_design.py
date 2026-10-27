from pathlib import Path
import json,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';s=(O/'metadata/GSE53276_samples.soft').read_bytes().decode('utf-8');out=[]
for b in s.split('^SAMPLE = ')[1:]:
 ls=b.splitlines();r={'accession':ls[0],'fields':[l for l in ls if l.startswith(('!Sample_title','!Sample_characteristics','!Sample_description','!Sample_data_processing','!Sample_platform_id'))]};table=b.split('!sample_table_begin\n',1);r['table_header']=table[1].splitlines()[0] if len(table)>1 else None;r['table_rows']=len(table[1].split('!sample_table_end')[0].splitlines())-1 if len(table)>1 else 0;out.append(r)
print('GSE53276',json.dumps(out,ensure_ascii=False)[:9500]);root=ET.fromstring((O/'literature/PMC4776110.xml').read_bytes());pars=[]
for e in root.iter('p'):
 t=' '.join(e.itertext())
 if any(q in t.lower() for q in ['gdes','cell lines','gse30601','expression signature']):pars.append(t)
(O/'metadata/INTRINSIC_DESIGN_54.json').write_text(json.dumps({'GSE53276':out,'DDR2_2016_methods':pars},ensure_ascii=False,indent=2),encoding='utf-8');print('DDR2_2016','\n'.join(pars)[:10500])
