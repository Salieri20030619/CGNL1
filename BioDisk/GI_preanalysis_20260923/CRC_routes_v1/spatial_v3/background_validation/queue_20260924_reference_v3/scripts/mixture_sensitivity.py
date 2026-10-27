import os,runpy,pathlib
os.environ['REFERENCE_FILE']='reference_normal_hep_only.npz'
os.environ['MIXTURE_OUTPUT']='V03_HEPREF_SENS'
runpy.run_path(str(pathlib.Path(__file__).with_name('mixture.py')),run_name='__main__')
