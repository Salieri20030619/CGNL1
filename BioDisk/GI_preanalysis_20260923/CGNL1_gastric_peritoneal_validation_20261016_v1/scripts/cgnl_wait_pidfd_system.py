import os,select,sys
try:
 fd=os.pidfd_open(int(sys.argv[1]));poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
except ProcessLookupError:pass
