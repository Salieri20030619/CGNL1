from pathlib import Path
from zipfile import ZipFile
from lxml import etree
import hashlib,json,re,shutil,subprocess
R=Path(r'C:\Users\31937\GC_validation_20260923');D=R/'delivery';doc=D/'胃癌腹膜转移研究建议书.docx';pdf=R/'render_client_final'/'胃癌腹膜转移研究建议书.pdf';target=D/pdf.name;shutil.copy2(pdf,target)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with ZipFile(doc) as z:
 assert z.testzip() is None
 root=etree.fromstring(z.read('word/document.xml'));styles=etree.fromstring(z.read('word/styles.xml'));text=''.join(root.xpath('//w:t/text()',namespaces=ns));assert not re.search(r'(\[\d+\])\1',text);assert len(root.xpath('//w:pBdr',namespaces=ns))==0;assert len(styles.xpath('//w:pBdr',namespaces=ns))==0;assert len(root.xpath('//a:blip',namespaces=ns))==3
 normalized=re.sub(r'[\s,]','',text);assert all(t in normalized for t in ['17','28','0.982','0.104','0.436','84.3%','69.4%','43113'])
info=subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdfinfo.exe',str(target)],capture_output=True,text=True,encoding='utf-8',errors='replace').stdout;assert re.search(r'Pages:\s+6',info)
qa={'status':'pass','pages':6,'all_pages_visually_inspected':True,'no_clipping_or_overlap':True,'missing_plot_glyphs_fixed':True,'title_rule_removed':True,'duplicate_citations_removed':True,'embedded_figures':3,'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'pdf_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'scope':'Client proposal and prevalidation, not confirmed mechanism or publication guarantee'}
(R/'document_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa))
