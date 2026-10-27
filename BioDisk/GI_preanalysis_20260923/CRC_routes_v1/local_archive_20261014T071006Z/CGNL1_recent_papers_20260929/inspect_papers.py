from pathlib import Path
import urllib.request,concurrent.futures,xml.etree.ElementTree as E,re
R=Path(r'C:/Users/31937/CGNL1_recent_papers_20260929')
ids=['PMC12267411','PMC12780015']
def f(pid):
 try:
  op=urllib.request.build_opener(urllib.request.ProxyHandler({}))
  with op.open('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pid+'/fullTextXML',timeout=25) as h:b=h.read()
  (R/(pid+'.xml')).write_bytes(b);root=E.fromstring(b)
  txt='\n'.join(''.join(x.itertext()) for x in root.findall('.//body//p'))
  (R/(pid+'.txt')).write_text(txt,encoding='utf-8')
  print('\nPAPER',pid,'SECTIONS', [''.join(t.itertext()) for t in root.findall('.//sec/title')])
  for x in root.findall('.//p'):
   s=''.join(x.itertext())
   if re.search(r'ethic|informed consent|generated|collected|cell culture|immunohist|western blot|publicly|cohorts|patient-derived|viability assay',s,re.I):print(s[:1400])
 except Exception as e:print(pid,'ERROR',str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(f,ids))
