from common import *
import subprocess,time,os
env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
for name,script in [('V04_DIAGNOSTIC','external_diagnostic.py'),('V04_BROAD_EP_SENS','external_broad_sensitivity.py')]:
 D=Q/name;D.mkdir(exist_ok=True)
 with (D/'run.log').open('w') as f:r=subprocess.run([sys.executable,str(Q/'scripts'/script),str(Q)],stdout=f,stderr=subprocess.STDOUT,env=env)
 dump(D/'EXIT.json',{'exit_code':r.returncode,'result_exists':(D/'RESULT.json').exists()})
 if r.returncode==0 and (D/'RESULT.json').exists():dump(D/'COMPLETE.json',{'exit_code':0})
 if r.returncode:raise RuntimeError(script+' failed')
