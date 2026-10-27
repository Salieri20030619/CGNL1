from pathlib import Path
import ctypes,os,json
b=Path(__file__).resolve().parents[1];m=b/'metadata';lib=ctypes.CDLL(None,use_errno=True);fd=lib.inotify_init1(0);assert fd>=0;wd=lib.inotify_add_watch(fd,os.fsencode(m),0x8|0x80);assert wd>=0
while not (m/'PILOT_DOWNLOAD_COMPLETE_451.json').exists():
 assert not (m/'LOCAL_RELAY_FAILED_451.json').exists(), 'Local relay failed; inspect failure marker'
 os.read(fd,65536)
os.close(fd);assert json.loads((m/'PILOT_DOWNLOAD_COMPLETE_451.json').read_text())['all_md5_verified'];print('RELAY_COMPLETE')
