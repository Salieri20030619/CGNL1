from pathlib import Path
import json,hashlib,datetime,tarfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];C=json.loads((R/'SPECIFICATION.json').read_text());D=R/'results'
assert json.loads((R/'PROCESS_EXIT.json').read_text())['exit_code']==0
res=pd.read_csv(D/'continuous_models_all_stages.tsv',sep='\t');mt=pd.read_csv(D/'within_patient_models_all_stages.tsv',sep='\t')
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
stages=[s['stage'] for s in C['raw_stages']];xx=np.arange(len(stages));colors=['#708090','#176b9b','#c86520']
fig,axes=plt.subplots(3,1,figsize=(12,10),sharex=True)
for model,color in zip(C['models'],colors):
 v=res[(res.branch=='raw')&(res.scope=='all')&(res.model==model)].set_index('stage').reindex(stages)
 axes[0].plot(xx,v.beta,'o-',label=model,color=color);axes[2].plot(xx,v.P,'o-',label=model,color=color)
for scope,color in [('all','#176b9b'),('external','#c86520')]:
 v=res[(res.branch=='raw')&(res.scope==scope)&(res.model=='M1_study')].set_index('stage').reindex(stages);axes[1].plot(xx,v.n_patients,'o-',label=scope,color=color)
axes[0].axhline(0,color='grey',lw=.8);axes[0].set_ylabel('Standardized association');axes[0].legend(ncol=3,fontsize=9)
axes[1].set_ylabel('Independent patient labels');axes[1].legend()
axes[2].set_ylabel('Raw P value');axes[2].set_yscale('log');axes[2].axhline(.05,color='firebrick',ls='--',label='P = 0.05');axes[2].set_xticks(xx,stages,rotation=35,ha='right');axes[2].legend()
fig.suptitle('CGNL1: broad-to-strict threshold ladder\nAll raw-count cohorts; all prespecified stages retained',fontsize=14);fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(R/'figures/continuous_threshold_ladder.png',dpi=170);fig.savefig(R/'figures/continuous_threshold_ladder.svg');plt.close(fig)
stages=[s['stage'] for s in C['match_stages']];xx=np.arange(len(stages));fig,axes=plt.subplots(3,1,figsize=(13,10),sharex=True)
for scope,color in [('all','#176b9b'),('external','#c86520')]:
 v=mt[mt.scope==scope].set_index('stage').reindex(stages);axes[0].plot(xx,v.effect,'o-',color=color,label=scope);axes[1].plot(xx,v.n_patients,'o-',color=color,label=scope);axes[2].plot(xx,v.P,'o-',color=color,label=scope)
axes[0].axhline(0,color='grey',lw=.8);axes[0].set_ylabel('Mean within-patient difference');axes[0].legend();axes[1].set_ylabel('Independent patient labels');axes[2].set_ylabel('Raw P value');axes[2].set_yscale('log');axes[2].axhline(.05,color='firebrick',ls='--');axes[2].set_xticks(xx,stages,rotation=45,ha='right')
fig.suptitle('Within-patient CGNL1 detection-associated junction state\nMatching, QC, identity and sample-size restrictions added step by step',fontsize=14);fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(R/'figures/within_patient_threshold_ladder.png',dpi=170);fig.savefig(R/'figures/within_patient_threshold_ladder.svg');plt.close(fig)
# Matching threshold crossings, including regained nominal significance.
trans=[]
for scope in mt.scope.unique():
 prev=None
 for x in mt[mt.scope==scope].to_dict('records'):
  if prev and np.isfinite(prev['P']) and np.isfinite(x['P']) and (prev['P']<.05)!=(x['P']<.05):trans.append(dict(scope=scope,from_stage=prev['stage'],to_stage=x['stage'],transition='gain' if x['P']<.05 else 'loss',P_before=prev['P'],P_after=x['P'],n_before=prev['n_patients'],n_after=x['n_patients']))
  prev=x
pd.DataFrame(trans,columns=['scope','from_stage','to_stage','transition','P_before','P_after','n_before','n_after']).to_csv(D/'within_patient_nominal_transitions.tsv',sep='\t',index=False)
files=[f for s in ['results','figures','scripts','logs'] for f in (R/s).rglob('*') if f.is_file()]+[R/'SPECIFICATION.json',R/'COMPLETE.json',R/'PROCESS_EXIT.json',R/'QA.json',R/'INTERPRETATION.json']
manifest=[dict(path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(R/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2));files.append(R/'DELIVERY_MANIFEST.json')
with tarfile.open(R/'threshold_ladder_delivery.tar.gz','w:gz') as t:
 for f in files:t.add(f,arcname=str(f.relative_to(R)))
print('DELIVERY_READY',len(files))
