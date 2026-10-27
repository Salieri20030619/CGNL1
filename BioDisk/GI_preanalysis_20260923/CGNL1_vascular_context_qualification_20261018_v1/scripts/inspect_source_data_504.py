from pathlib import Path
import zipfile,xml.etree.ElementTree as ET,json,hashlib,re
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
p=B/'raw/source_linkage_503/source_data.xlsx'
assert hashlib.md5(p.read_bytes()).hexdigest()=='d97e1b76f15ae20dea087e1650ff4c7f'
n={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(p) as z:
    shared=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        shared=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',n)]
    rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    report=[]
    for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',n):
        target=rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        path=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        root=ET.fromstring(z.read(path));texts=[];formulas=0;cells=0
        for c in root.findall('.//s:sheetData/s:row/s:c',n):
            cells+=1
            formulas+=c.find('s:f',n) is not None
            v=c.find('s:v',n)
            if c.attrib.get('t')=='s': value=shared[int(v.text)]
            elif c.attrib.get('t')=='inlineStr': value=''.join(c.find('s:is',n).itertext())
            else: continue
            if value.strip(): texts.append(dict(cell=c.attrib['r'],text=value))
        report.append(dict(sheet=sheet.attrib['name'],populated_cells=cells,formula_cells=formulas,text_cells=texts))
    out=dict(source=str(p),source_md5=hashlib.md5(p.read_bytes()).hexdigest(),sheets=report,read_only=True)
    (B/'metadata/SOURCE_DATA_INVENTORY_504.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    for s in report:
        print(json.dumps(dict(sheet=s['sheet'],headings=s['text_cells'][:3],RNA_identifier_hits=[x for x in s['text_cells'] if re.search('donor|KR[-_]|GSM|RNA|sample',x['text'],re.I)])))

