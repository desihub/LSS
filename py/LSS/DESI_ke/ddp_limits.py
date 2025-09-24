import os
import sys
import argparse
import numpy             as np
import matplotlib.pyplot as plt
import cosmo             as cosmo
import astropy.io.fits   as fits

from   astropy.table     import Table
from   smith_kcorr       import GAMA_KCorrection, GAMA_KCorrection_color
from   rest_gmr          import smith_rest_gmr
from   tmr_ecorr         import tmr_ecorr, tmr_q
from   abs_mag           import abs_mag
from   data.ke_params    import *
from   findfile          import findfile, fetch_header
from   config            import Configuration


parser   = argparse.ArgumentParser(description='Gen kE DDP limit curves')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
parser.add_argument('-s', '--survey', help='Select survey', default='gama')

args     = parser.parse_args()
log      = args.log
survey   = args.survey.lower()

config = Configuration(args.config)
config.update_attributes('ddp_limits', args)
config.write()
    
if log:
    logfile = findfile(ftype='ddp_limit', dryrun=False, survey=survey, log=True)

    print(f'Logging to {logfile}')

    sys.stdout = open(logfile, 'w')

kcorr_r  = GAMA_KCorrection(band='R')
kcorr_RG = GAMA_KCorrection_color()

# To be looped over for a total of 28 = 7 (rest color) x 2 (Qall, Qcolor) x 2 (magnitude) curves.
gmrs_0p1 = np.array([0.131, 0.298, 0.443, 0.603, 0.785, 0.933, 1.067])  
gmrs_0p0 = np.array([0.158, 0.298, 0.419, 0.553, 0.708, 0.796, 0.960])

# bright and faint limits.   
rlims    = [fetch_header(ftype='gold', name='RMAX', survey=survey),\
            fetch_header(ftype='gold', name='RLIM', survey=survey)]

root     = os.environ['GOLD_DIR'] + f'/ddrp_limits/'

if not os.path.isdir(root):
    print('Creating {}'.format(root))

    os.makedirs(root)

count    = 0

zs = mus = None

for rlim in rlims:
    print('----------------------------------')

    rs = rlim * np.ones_like(zs)

    for aall, all_type in zip([True, False], ['QALL', 'QCOLOR']):
        for gmr_0P1 in gmrs_0p1:
            opath    = findfile(ftype='ddp_limit', dryrun=False, survey=survey, ddp_count=count)

            if args.nooverwrite & os.path.isfile(opath):
                hdul = fits.open(opath)
                hdr  = hdul[0].header  # the primary HDU header

                assert 'SURVEY' in hdr

                if hdr['SURVEY'].upper() == survey.upper():
                    print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
                
                    count += 1

                    continue

            if (zs is None) | (mus is None):
                zs   = np.arange(0.01, 0.6, 0.01)
                mus  = cosmo.distmod(zs)

            gmr_0P1  = gmr_0P1 * np.ones_like(zs)
            gmr_0P0  = kcorr_RG.rest_gmr_nonnative(gmr_0P1)

            ks       = kcorr_r.k_nonnative_zref(0.0, zs, gmr_0P1)
            es       = tmr_ecorr(zs, gmr_0P0, aall=aall)
            Mrs_0P0  = abs_mag(rs, mus, ks, es)

            dat      = Table(np.c_[zs, ks, es, Mrs_0P0], names=['Z', 'K', 'E', 'M0P0_{}'.format(all_type)])
            dat.meta = {'RLIM': rlim, 'ALL': aall, 'GMR_0P1': gmr_0P1[0], 'GMR_0P0': gmr_0P0[0], 'SURVEY': survey}
            
            dat.write(opath, format='fits', overwrite=True)
            
            count   += 1

            print()
            print('Solved for {} {} {}: {}'.format(rlim, all_type, gmr_0P1[0], opath))

print('Done.')

if log:
    sys.stdout.close()
