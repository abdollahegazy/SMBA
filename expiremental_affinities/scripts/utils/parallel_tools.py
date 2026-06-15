import os
import signal

def _worker_init():
    os.setpgrp()
    signal.signal(signal.SIGINT, signal.SIG_IGN)
