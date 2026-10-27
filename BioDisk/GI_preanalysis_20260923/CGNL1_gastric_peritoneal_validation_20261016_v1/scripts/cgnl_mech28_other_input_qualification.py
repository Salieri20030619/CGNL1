from pathlib import Path
import gzip,json,zipfile,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';out={}
for name in ['GSE121372','GSE120303']:
 p=O/'metadata'/(name+'_samples.soft');s=p.read_bytes().decode('utf-8');samples=[]
 for b in s.split('^SAMPLE = ')[1:]:
  ls=b.splitlines();keep=[l for l in ls if l.startswith(('!Sample_title =','!Sample_characteristics_ch1 =','!Sample_platform_id =','!Sample_description =','!Sample_data_processing =','!Sample_treatment_protocol_ch1 ='))];samples.append({'accession':ls[0],'fields':keep})
 out[name]=samples;print(name,json.dumps(samples,ensure_ascii=False)[:6500])
for p in (O/'raw').glob('GSE121372*.gz'):
 try:
  b=gzip.decompress(p.read_bytes());t=b.decode('utf-8');out[p.name]={'crc':True,'header':t.splitlines()[:3]};print(p.name,t[:800])
 except Exception as e:out[p.name]={'crc':False,'error':str(e)};print(p.name,str(e))
for p in (O/'metadata').glob('41467_2023_36310*.xlsx'):
 with zipfile.ZipFile(p) as z:
  hits=[]
  for n in z.namelist():
   if n.endswith('.xml'):
    text=z.read(n).decode('utf-8')
    if 'CGNL1' in text or 'ENSG00000128849' in text:hits.append(n)
  out[p.name]={'target_identifier_xml_hits':hits,'full_matrix_present':False,'reason':'small figure source-data and signature workbooks, not191987cell transcriptome matrix'}
(O/'metadata/INPUT_QUALIFICATION_28.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
