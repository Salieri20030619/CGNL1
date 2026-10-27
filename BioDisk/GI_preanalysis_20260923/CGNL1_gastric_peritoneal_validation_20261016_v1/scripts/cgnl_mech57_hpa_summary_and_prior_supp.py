from pathlib import Path
import json,xml.etree.ElementTree as ET,collections,requests
O=Path.cwd()/'mechanism_investigation_v2';root=ET.fromstring((O/'metadata/HPA_CGNL1.xml').read_bytes());out=[]
for ab in root.findall('.//antibody'):
 for te in ab.findall('tissueExpression'):
  if te.attrib.get('assayType')!='cancer':continue
  for d in te.findall('data'):
   if d.findtext('tissue')!='Stomach cancer':continue
   pats=[]
   for p in d.findall('patient'):pats.append({'patient_id':p.findtext('patientId'),'staining':p.findtext("level[@type='staining']"),'intensity':p.findtext("level[@type='intensity']"),'location':p.findtext('location'),'image_urls':[e.text for e in p.findall('.//imageUrl')]})
   out.append({'antibody':ab.attrib['id'],'n_patients':len(pats),'counts':dict(collections.Counter(p['staining'] for p in pats)),'patients':pats,'antibody_validation':[ET.tostring(e,encoding='unicode') for e in ab if 'validation' in e.tag.lower()]})
(O/'results/HPA_GASTRIC_PROTEIN_SUMMARY.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('HPA',json.dumps([{k:v for k,v in r.items() if k not in ['patients','antibody_validation']} for r in out]))
root=ET.fromstring((O/'literature/PMC4776110.xml').read_bytes());links=[e.attrib for e in root.iter() if e.tag in ['media','supplementary-material','ext-link'] and any(x in str(e.attrib) for x in ['xls','supp','srep22371'])];print('PRIOR_SUPP',links)
u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4776110/supplementaryFiles'
try:
 a=requests.get(u,timeout=(15,90));a.raise_for_status();(O/'metadata/DDR2_2016_supplement.zip').write_bytes(a.content);print('SUPP_DOWNLOAD',len(a.content))
except Exception as e:print('SUPP_ERROR',repr(e))
