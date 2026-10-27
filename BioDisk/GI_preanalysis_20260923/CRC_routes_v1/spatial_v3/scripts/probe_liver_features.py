import pathlib,sys,tarfile,gzip,pandas as pd,json
O=pathlib.Path(sys.argv[1]);plan=json.loads((O/'SCREEN_PLAN.json').read_text())
for name in ['IX_56689','M-ST-11']:
 with tarfile.open(O/'raw/liver'/(name+'.tar.gz')) as t:
  f=[m for m in t.getnames() if m.endswith('features.tsv.gz')][0];x=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(f)),sep='\t',header=None);sy=set(x[1].astype(str));print(name,len(x));print('FOCUSED',[(g,g in sy) for g in plan['candidates']+['KRT8','KRT18','KRT19','KRT20','EPCAM','CDH1']]);print('PROGRAMS',[(n,sorted(sy&set(gs))) for n,gs in plan['programs'].items()])
