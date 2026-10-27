import pathlib,xml.etree.ElementTree as ET,re
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2'
root=ET.parse(B/'01_metadata'/'GC_228598.xml').getroot()
for e in root.iter('table-wrap'):
 s=' '.join(e.itertext());print(s[:6000])
for e in root.iter('supplementary-material'):print(ET.tostring(e,encoding='unicode')[:2000])
# Read-only inspect of public interface configuration, no authentication or mutations.
s=(B/'01_metadata'/'iprox_project.html').read_text()
for line in s.splitlines():
 if 'jsonServer' in line or 'config' in line.lower() or ('script' in line and 'src' in line):print(line.strip()[:400])
