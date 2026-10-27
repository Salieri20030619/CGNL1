import pathlib,sys,tarfile,json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
O=pathlib.Path(sys.argv[1]);R=O/'results/screen';a=pd.read_csv(R/'spot_target_coordinates.tsv.gz',sep='\t');jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text());pal={'Hepatic lobule':'#6BAED6','Liver macrometastasis tumor':'#C44E52','Liver macrometastasis stroma':'#999999','Liver micrometastasis tumor':'#AD3BA4','Liver micrometastasis stroma':'#CCA2CD'}
for job in jobs:
 name=job['sample'];v=a[a['sample']==name];xy=v[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();macro=(v.Layer3=='Liver macrometastasis tumor').to_numpy();hep=(v.Layer3=='Hepatic lobule').to_numpy();good=v.QC_pass.to_numpy();d=v.distance_to_hepatic_um.to_numpy()
 with tarfile.open(O/'raw/liver'/(name+'.tar.gz')) as t:
  ns=t.getnames();sf=json.load(t.extractfile([p for p in ns if p.endswith('scalefactors_json.json')][0]));im=Image.open(t.extractfile([p for p in ns if p.endswith('tissue_hires_image.png')][0])).copy()
 fig,axs=plt.subplots(1,3,figsize=(14,4.8));axs[0].imshow(im);scale=sf['tissue_hires_scalef'];axs[0].scatter(xy[:,0]*scale,xy[:,1]*scale,c=[pal.get(z,'#DBCE9A') for z in v.Layer3],s=1,alpha=.25);axs[0].set_title('Pathology overlay',fontsize=10)
 sc=axs[1].scatter(xy[:,0],xy[:,1],c=np.where(good,v.CGNL1_log2CPM,np.nan),s=2,cmap='viridis',vmin=0,vmax=12);axs[1].invert_yaxis();axs[1].set_title('CGNL1 log2(CPM+1), QC spots',fontsize=10);fig.colorbar(sc,ax=axs[1],shrink=.65)
 color=np.full(len(v),'#E5E5E5',dtype=object);color[hep]='#A6CEE3';color[macro&good&(d<=300+1e-6)]='#E66101';color[macro&good&(d>=600-1e-6)]='#5E3C99';axs[2].scatter(xy[:,0],xy[:,1],c=color,s=2);axs[2].invert_yaxis();axs[2].set_title('Tumor: near (orange), far (purple)',fontsize=10)
 for ax in axs:ax.set_aspect('equal');ax.axis('off')
 fig.suptitle(f'Patient {job["patient"]} | {name} | '+('Recurrence' if v.recurrence.iloc[0] else 'Non-recurrence'),fontsize=13);fig.subplots_adjust(left=.02,right=.98,top=.86,bottom=.10,wspace=.08);fig.text(.5,.025,'Pathology: tumor red / hepatic blue / stroma grey.   Near <=300 um; far >=600 um from hepatic spot centers.',ha='center',fontsize=9);fig.savefig(R/'figures'/(name+'.png'),dpi=140);plt.close(fig)
files=[R/'figures'/(j['sample']+'.png') for j in jobs]
for start in range(0,len(files),4):
 im=Image.new('RGB',(1960,720),(255,255,255))
 for k,p in enumerate(files[start:start+4]):
  x=Image.open(p).convert('RGB');x.thumbnail((980,345));im.paste(x,((k%2)*980,(k//2)*360))
 im.save(R/'figures'/f'contact_{start//4+1}.jpg',quality=92)
print('24 spatial maps regenerated; statistics unchanged',flush=True)
