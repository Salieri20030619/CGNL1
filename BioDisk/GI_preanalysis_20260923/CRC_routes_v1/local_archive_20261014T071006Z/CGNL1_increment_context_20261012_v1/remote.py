import base64, pathlib, subprocess, sys
src=pathlib.Path(sys.argv[1]).read_bytes()
cmd='echo '+base64.b64encode(src).decode()+' | base64 -d | '+(sys.argv[2] if len(sys.argv)>2 else 'python3')
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10',cmd],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
sys.stdout.buffer.write(p.stdout);sys.stderr.buffer.write(p.stderr);sys.exit(p.returncode)
