import importlib,os,psutil,json
for n in ['scanpy','anndata','sklearn','igraph','leidenalg','harmonypy']:
 try:m=importlib.import_module(n);print(n,getattr(m,'__version__','available'))
 except Exception as e:print(n,'UNAVAILABLE',str(e))
print('RAM_GB',psutil.virtual_memory().available/1e9)
