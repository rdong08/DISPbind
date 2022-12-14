import sys
import os
import os.path
import math
import time
from collections import defaultdict
from functools import wraps
import pysam
try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans


def which(program):
    '''
    Check the path of external programs
    '''
    def is_executable(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    for path in os.environ["PATH"].split(os.pathsep):
        progpath = os.path.join(path, program)
        if is_executable(progpath):
            return progpath
    return None


def logger(fn):
    '''
    Record parameters and runming time
    '''
    @wraps(fn)
    def wrapper(*args, **kwargs):
        print(kwargs['command'])
        local_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        print('Start DISPbind %s at %s' % (kwargs['name'], local_time))
        fn(*args)
        local_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        print('End DISPbind %s at %s' % (kwargs['name'], local_time))
    return wrapper


