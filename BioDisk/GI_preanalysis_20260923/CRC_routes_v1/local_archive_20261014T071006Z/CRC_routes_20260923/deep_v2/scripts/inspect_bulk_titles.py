import gzip,pathlib,re
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2')
s=gzip.open(O/'metadata/GSE50760.soft.gz','rt',encoding='utf-8').read()
for c in s.split('^SAMPLE = ')[1:]:print(c.splitlines()[0],re.search(r'!Sample_title = (.+)',c).group(1),re.findall(r'!Sample_characteristics_ch1 = tissue: (.+)',c))
