import pathlib,os,signal,json
p=pathlib.Path('/proc/3143461/cmdline')
if p.exists():
 b=p.read_bytes()
 expected=b'/spatial_v3/background_validation/final_exploration_20260924/scripts/final_source_state.py'
 if expected not in b:raise RuntimeError('PID does not match superseded project task')
 os.kill(3143461,signal.SIGTERM);print('Stopped only superseded source-state worker 3143461')
else:print('Previous worker already ended')
