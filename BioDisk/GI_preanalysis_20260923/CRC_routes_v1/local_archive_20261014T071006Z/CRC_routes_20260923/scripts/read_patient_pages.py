import pathlib,pypdf,json,collections,subprocess
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');r=pypdf.PdfReader(O/'metadata/sciadv.adf5464_sm.pdf')
for i in [10,11,12]:print('PAGE',i+1,r.pages[i].extract_text(extraction_mode='layout'))
