import pathlib,h5py
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for name in ['GSE294559_TP001.h5','GSE294559_TO001.h5']:
 with h5py.File(B/'02_raw'/name,'r') as f:
  m=f['matrix'];print(name,list(m.keys()),list(m['features'].keys()));print('SHAPE',m['shape'][:]);print('FEATURE',m['features']['name'][:5],m['features']['feature_type'][:5],m['features']['name'][-3:])
