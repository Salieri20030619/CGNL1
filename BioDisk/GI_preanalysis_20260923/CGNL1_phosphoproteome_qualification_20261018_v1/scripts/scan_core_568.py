from pathlib import Path
import zipfile,xml.etree.ElementTree as E,json,re,datetime,posixpath
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');out=B/'results/core_qualification_568';out.mkdir(exist_ok=True)
files=[('HCC_protein',B/'raw/intermittent_564/HCC_verified/mmc3.xlsx'),('HCC_phospho',B/'raw/intermittent_564/HCC_verified/mmc6.xlsx'),('CRC',B/'raw/intermittent_564/CRC_verified/MCO2-6-e70492-s011.xlsx'),('CCA',B/'raw/intermittent_564/CCA_static.verified')]
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};rel='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
for tag,path in files:
 rec={'source':str(path),'sheets':[]}
 with zipfile.ZipFile(path) as z:
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   with z.open('xl/sharedStrings.xml') as f:
    for event,el in E.iterparse(f,events=('end',)):
     if el.tag.endswith('}si'):strings.append(''.join(x.text or '' for x in el.iter() if x.tag.endswith('}t')));el.clear()
  links={x.get('Id'):x.get('Target') for x in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))};sheets=E.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',ns)
  def value(c):
   v=c.find('m:v',ns);t=c.get('t')
   if t=='inlineStr':return ''.join(x.text or '' for x in c.iter() if x.tag.endswith('}t'))
   if v is None:return None
   if t=='s':return strings[int(v.text)]
   if t in ['str','e','b']:return v.text
   try:return float(v.text)
   except (ValueError,TypeError):return v.text
  for sheet in sheets:
   name=sheet.get('name');target=links[sheet.get(rel)];member=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);heads=[];hits=[];n=0
   with z.open(member) as f:
    for event,el in E.iterparse(f,events=('end',)):
     if not el.tag.endswith('}row'):continue
     n+=1;cells=el.findall('m:c',ns);first=[value(c) for c in cells[:5]];matched=any(isinstance(x,str) and ('CGNL1' in x.upper() or 'Q0VF96' in x.upper()) for x in first)
     if n<=3 or matched or (len(cells)<=2 and n<=15):
      row={'row':int(el.get('r',n)),'cells':{c.get('r'):value(c) for c in cells if value(c) is not None}}
      if n<=3 or (len(cells)<=2 and n<=15):heads.append(row)
      if matched:hits.append(row)
     el.clear()
   rec['sheets'].append({'name':name,'rows':n,'headers':heads,'target_rows':hits});print(tag,name,'rows',n,'target_hits',len(hits),flush=True)
 rec['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(out/(tag+'.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8')
