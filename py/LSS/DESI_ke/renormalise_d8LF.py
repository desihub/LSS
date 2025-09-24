import numpy as np

from astropy.table import Table
from ddp import tmr_DDP1


def renormalise_d8LF(idx, cat, fdelta, fdelta_zeropoint, self_count=False):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''

    print('Renormalising (IVMAX) LF with log10(fdelta)={:.6e} for tier {}'.format(np.log10(fdelta), idx))
    
    cat = Table(cat, copy=True)

    cols = ['PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR']
    
    if 'PHI_IVMAX_JK' in cat.dtype.names:
        print('Including JK estimates in renormalisation.')

        cols += ['PHI_IVMAX_JK', 'PHI_IVMAX_ERROR_JK']

    for col in cols:
        cat[col] /= fdelta
    
    if self_count:
        print('Applying log10|self-count correction| of {:.6f}'.format(np.log10(fdelta / fdelta_zeropoint)))
            
        # tmr_DDP1: [-21.8, -20.1] 
        is_ddp1 = (cat['MEDIAN_M'] > tmr_DDP1[0]) & (cat['MEDIAN_M'] < tmr_DDP1[1])

        for col in cols:
            cat[col][is_ddp1] *= (fdelta / fdelta_zeropoint)

    return  cat
