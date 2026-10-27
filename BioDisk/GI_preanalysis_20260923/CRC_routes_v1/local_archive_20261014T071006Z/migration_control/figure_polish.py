from pathlib import Path
import subprocess
p=Path(r'D:\CGNL1\_transfer\05_report.py');s=p.read_text(encoding='utf-8')
s=s.replace('import matplotlib.pyplot as plt','import matplotlib.pyplot as plt\nfrom matplotlib.ticker import MaxNLocator')
s=s.replace("axes[0].plot([0,1],y,color='#327f93',alpha=.7,marker='o',ms=3)","axes[0].plot([0,1],y,color='#24778d' if g.tissueSite.iloc[0]=='Liver' else '#cf8646',alpha=.7,marker='o',ms=3)")
s=s.replace("xytext=(3,3),textcoords='offset points',fontsize=7", "xytext={'P03':(-17,8),'P30':(5,-12),'P13':(-15,-12),'P28':(4,8)}.get(row.donor_id,(3,3)),textcoords='offset points',fontsize=7")
s=s.replace("fig.suptitle('SCP2079", "for ax in axes[1:]:ax.xaxis.set_major_locator(MaxNLocator(5))\nfig.suptitle('SCP2079")
compile(s,str(p),'exec');p.write_text(s,encoding='utf-8')
flags=subprocess.CREATE_NO_WINDOW
subprocess.run(['scp','-o','BatchMode=yes',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/SCP2079_paired_20261014_v1/scripts/'],check=True,capture_output=True,creationflags=flags)
code='''
from pathlib import Path
import subprocess,json,datetime
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/SCP2079_paired_20261014_v1')
with (R/'logs/05_report.py.log').open('w',encoding='utf-8') as f:p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(R/'scripts/05_report.py')],stdout=f,stderr=subprocess.STDOUT)
r={'script':'05_report.py','exit_code':p.returncode,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'logs/05_report.py.exit.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r));assert p.returncode==0
'''
r=subprocess.run(['ssh','-o','BatchMode=yes','sgs10','python3 -'],input=code.encode('utf-8'),capture_output=True,creationflags=flags);assert r.returncode==0,r.stderr.decode('utf-8');print(r.stdout.decode('utf-8'))