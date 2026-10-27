import pathlib,sys,importlib.util
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');sys.path.insert(0,str(O/'vendor'))
for x in ['pybiomart','pyreadr','pytoml','session_info','infercnvpy']:print(x,bool(importlib.util.find_spec(x)))
try:import infercnvpy;print('IMPORT_OK',infercnvpy.__version__)
except Exception as e:print(type(e).__name__,str(e))
