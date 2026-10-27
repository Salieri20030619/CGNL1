import pathlib,gzip,json,re,sys
sys.stdout.reconfigure(encoding='utf-8');O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\metadata')
j=json.loads((O/'zenodo19043057.json').read_text());print('ZENODO',j.get('links'));print([(f['key'],f['size'],f['links']) for f in j.get('files',[])])
for fn in O.glob('*.soft.gz'):
 s=gzip.open(fn,'rt',encoding='utf-8').read();print('\n'+fn.name)
 print('\n'.join(x for x in s.splitlines() if x.startswith(('!Series_title','!Series_pubmed','!Series_overall_design','!Series_supplementary_file'))))
 for chunk in s.split('^SAMPLE = ')[1:]:
  lines=chunk.splitlines();print(lines[0],'; '.join(x.split(' = ',1)[-1] for x in lines if x.startswith(('!Sample_title','!Sample_characteristics_ch1','!Sample_supplementary_file'))))
