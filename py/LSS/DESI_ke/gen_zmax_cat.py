import os
import sys
import time
import argparse
import runtime
import numpy as np

from   cosmo import distmod, volcom
from   smith_kcorr import GAMA_KCorrection
from   tmr_ecorr import tmr_ecorr
from   scipy.optimize import brentq, minimize
from   astropy.table import Table
from   functools import partial
from   multiprocessing import Pool
from   findfile import findfile, overwrite_check, write_desitable, fetch_header
from   config import Configuration

kcorr_r = GAMA_KCorrection(band='R')

def theta(z, rest_gmr_0p1, rest_gmr_0p0, aall=False):
    z            = np.atleast_1d(z)
    rest_gmr_0p1 = np.atleast_1d(rest_gmr_0p1)
    rest_gmr_0p0 = np.atleast_1d(rest_gmr_0p0)
    
    result       = distmod(z) + kcorr_r.k_nonnative_zref(0.0, z, rest_gmr_0p1) + tmr_ecorr(z, rest_gmr_0p0, aall=aall)

    return  result[0]
    
def solve_theta(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall=False):
     def diff(x):
          return theta(x, rest_gmr_0p1, rest_gmr_0p0, aall=aall) - thetaz - dr

     def absdiff(x):
        return np.abs(diff(x))

     warn = 0

     try:
        result = brentq(diff, 1.e-3, 1.6)

     except ValueError as VE:
        warn = 1

        # Brent method fails, requires sign change across boundaries.                                                                                          
        result = minimize(absdiff, 1.)

        if result.success:
            result = result.x[0]

        else:
             warn = 2
             result = -99.

     return  result, warn

def zmax(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs, aall=False, debug=True):
   result = []
   start  = time.time()

   if debug:
        print('Solving for zlimit.')

   '''
   for i, (rest_gmr_0p1, rest_gmr_0p0, theta_z, dr) in enumerate(zip(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs)):
        interim, warn = solve_theta(rest_gmr_0p1, rest_gmr_0p0, theta_z, dr, aall=aall)

        result.append([interim, warn])

        if (i % 500 == 0) & debug:
             runtime = (time.time() - start) / 60.

             print('{:.3f}% complete after {:.2f} mins.'.format(100. * i / len(theta_zs), runtime))
   '''
   with Pool(processes=14) as pool:
       arglist = list(zip(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs))
       result  = pool.starmap(partial(solve_theta, aall=aall), arglist)
   
   result = np.array(result)

   return  result[:,0], result[:,1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gen zmax cat.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')  
    parser.add_argument('-a', '--aall',   help='All Q, no red/blue split.', action='store_true')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--theta_def',    help='Specifier for definition of theta', default='Z_THETA_QCOLOR')
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    
    args      = parser.parse_args()
    log       = args.log
    aall      = args.aall
    dryrun    = args.dryrun
    survey    = args.survey.lower()
    theta_def = args.theta_def

    config = Configuration(args.config)
    config.update_attributes('zmax', args)
    config.write()

    rlim      = fetch_header(ftype='gold', name='RLIM', survey=survey)
    rmax      = fetch_header(ftype='gold', name='RMAX', survey=survey)

    start     = time.time()

    if log:
        logfile = findfile(ftype='zmax', dryrun=False, survey=survey, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')

    print('Assuming {:.4f} < r < {:.4f}'.format(rmax, rlim))
    print('Assuming Q ALL = {}'.format(aall))
    
    fpath  = findfile(ftype='kE',   dryrun=dryrun, survey=survey)
    opath  = findfile(ftype='zmax', dryrun=dryrun, survey=survey)
    
    if args.nooverwrite:
        overwrite_check(opath)

    print('Reading {}.'.format(fpath))
        
    dat = Table.read(fpath)
    dat.pprint()

    dat['DELTA_DETMAG_FAINT'] = rlim - dat['DETMAG']

    print('Solving for {} bounding curve'.format(rlim))
    
    zmaxs, warn = zmax(dat['REST_GMR_0P1'],\
                       dat['REST_GMR_0P0'],\
                       dat[theta_def],\
                       dat['DELTA_DETMAG_FAINT'],\
                       aall=aall,\
                       debug=True)

    dat['ZMAX']      = zmaxs
    dat['ZMAX_WARN'] = warn

    dat['DELTA_DETMAG_BRIGHT'] = rmax - dat['DETMAG']

    print('Solving for {} bounding curve'.format(rmax))
    
    zmins, warn = zmax(dat['REST_GMR_0P1'],\
                       dat['REST_GMR_0P0'],\
                       dat[theta_def],\
                       dat['DELTA_DETMAG_BRIGHT'],\
                       aall=aall,\
                       debug=True)

    dat['ZMIN']      = zmins
    dat['ZMIN_WARN'] = warn

    dat.meta['THETA_DEF'] = theta_def

    dat['VMAX']  = volcom(dat['ZMAX'], dat.meta['AREA'])
    dat['VMAX'] -= volcom(dat['ZMIN'], dat.meta['AREA'])

    print('Writing {}.'.format(opath))

    dat.pprint()

    write_desitable(opath, dat)

    nwarn   = (dat['ZMAX_WARN'].data > 0) | (dat['ZMIN_WARN'].data > 0)
    nwarn   = np.count_nonzero(nwarn)

    print(f'WARNING:  zmax/min warnings triggered on {nwarn} galaxies.')

    runtime = (time.time() - start) / 60.

    print('\n\nDone in {} mins.\n\n'.format(runtime))

    if log:
        sys.stdout.close()
