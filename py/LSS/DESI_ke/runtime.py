import os
import sys
import time
import warnings

from   astropy.io.fits.verify import VerifyWarning
from   astropy.utils.metadata import MergeConflictWarning

# Suppress verify warnings, e.g. HIERARCH card length. 
warnings.simplefilter('ignore', category=VerifyWarning)
warnings.simplefilter('ignore', category=MergeConflictWarning)

def sizeofMB(xx):
    return  sys.getsizeof(xx)

def calc_runtime(start, log=None, memuse=False, xx=None):
    runtime  = time.time() - start
    runtime /= 60.

    if memuse:
        # COSMA5 has 302 compute nodes, each with 128GB RAM and 16 cores
        # https://www.dur.ac.uk/icc/cosma/cosma5/
            
        import psutil
        process = psutil.Process(os.getpid())
        memuse  = process.memory_info().rss / 1.e6 # [MB]
        memuse  = ' (memuse {}: {:.2f}MB)'.format(os.getpid(), memuse)

    else:
        memuse  = ''

    if log != None:
        msg = '{} after {:.4f} mins{}.'.format(log, runtime, memuse)
        
        if xx is not None:
            msg = msg.replace('Writing', 'Writing ({}MB)'.format(sizeofMB(xx)))
            msg = msg.replace('Reading', 'Reading ({}MB)'.format(sizeofMB(xx)))

        print(msg)
        
    return  runtime


if __name__ == '__main__':
    start = time.time()

    runtime  = calc_runtime(start, 'Read randoms')
