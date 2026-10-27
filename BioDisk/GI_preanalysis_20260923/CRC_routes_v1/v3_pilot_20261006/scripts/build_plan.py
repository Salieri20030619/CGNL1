from pathlib import Path
import json,hashlib,datetime
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from xml.sax.saxutils import escape
b=Path(__file__).resolve().parent
o=b/'output/pdf';o.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('Chinese',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0))
body=ParagraphStyle('body',fontName='Chinese',fontSize=10.5,leading=17,spaceAfter=12,wordWrap='CJK',textColor=HexColor('#223547'))
title=ParagraphStyle('title',parent=body,fontSize=20,leading=28,spaceAfter=24,textColor=HexColor('#084b61'))
small=ParagraphStyle('small',parent=body,fontSize=8,leading=12)
sections=json.loads((b/'plan_sections.json').read_text(encoding='utf-8'))
def footer(c,doc):
 c.setStrokeColor(HexColor('#cad7df'));c.line(42,44,A4[0]-42,44)
 c.setFont('Chinese',8);c.setFillColor(HexColor('#526b7d'))
 c.drawString(42,30,'CGNL1 | v3.0 | 2026-10-06 | 分析前冻结')
 c.drawRightString(A4[0]-42,30,str(doc.page))
story=[]
for i,(head,paras) in enumerate(sections):
 if i:story.append(PageBreak())
 story.append(Paragraph(escape(head),title))
 for s in paras:
  style=small if s.startswith('[') else body
  story.append(Paragraph(escape(s).replace('\n','<br/>'),style))
pdf=o/'CGNL1_research_plan_v3.0.pdf'
SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=47,bottomMargin=60).build(story,onFirstPage=footer,onLaterPages=footer)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
files=['plan_v3.json','plan_sections.json','output/pdf/CGNL1_research_plan_v3.0.pdf']
freeze={'version':'3.0','frozen_at_utc':now,'new_pilot_target_results_seen':False,'historical_results_known':True,'files':{f:hashlib.sha256((b/f).read_bytes()).hexdigest() for f in files}}
(b/'PLAN_FREEZE.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(freeze,ensure_ascii=False))
