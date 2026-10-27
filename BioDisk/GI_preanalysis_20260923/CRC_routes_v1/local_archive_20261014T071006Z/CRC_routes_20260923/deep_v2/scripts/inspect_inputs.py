import pathlib,gzip,re,json
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');s=gzip.open(O/'metadata/GSE50760.soft.gz','rt',encoding='utf-8').read();print('\n'.join(x for x in s.splitlines() if x.startswith(('!Series_supplementary','!Series_pubmed'))))
for c in s.split('^SAMPLE = ')[1:4]:print('\n'.join(x for x in c.splitlines() if x.startswith(('!Sample_title','!Sample_characteristics','!Sample_supplementary','!Sample_data_processing'))))
j=json.loads((O/'metadata/infercnvpy_pypi.json').read_text());print('INFER',j['info']['version'],j['info']['requires_dist']);print([(x['filename'],x['url']) for x in j['urls']])
print((O/'metadata/GSE50760.files.html').read_text())
