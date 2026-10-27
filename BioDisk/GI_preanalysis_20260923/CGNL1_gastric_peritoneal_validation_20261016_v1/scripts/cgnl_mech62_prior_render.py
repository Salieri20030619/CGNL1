from pathlib import Path
import zipfile,subprocess,json
O=Path.cwd()/'mechanism_investigation_v2'
with zipfile.ZipFile(O/'metadata/DDR2_2016_supplement.zip') as z:(O/'literature/srep22371-f1.jpg').write_bytes(z.read('srep22371-f1.jpg'))
subprocess.run(['pdftoppm','-f','5','-l','5','-r','120','-png','-singlefile',str(O/'literature/srep22371-s1.pdf'),str(O/'literature/DDR2_supp_page5')],check=True)
print('Rendered originalFigure1 andsupplementpage5 forvisualverification')
