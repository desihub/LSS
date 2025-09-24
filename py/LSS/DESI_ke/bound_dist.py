import os
import gc
import sys
import tqdm
import time
import fitsio
import argparse
import multiprocessing
import numpy             as np
import astropy.io.fits   as fits
import matplotlib.pyplot as plt

from   scipy.spatial   import KDTree
from   astropy.table   import Table
from   multiprocessing import Pool
from   runtime         import calc_runtime
from   findfile        import findfile, overwrite_check, call_signature
from   bitmask         import lumfn_mask, consv_mask
from   config          import Configuration

'''
Script to calculate the maximum distance [Mpc/h] of each random from the boundary. 
'''

np.random.seed(314)

parser = argparse.ArgumentParser(description='Find boundary distance for all randoms in a specified field..')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey.', default='gama')
parser.add_argument('--prefix', help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
parser.add_argument('--nproc', type=int, help='Number of processors', default=12)
parser.add_argument('--realz', type=int, help='Realisation', default=0)

args   = parser.parse_args()
log    = args.log
field  = args.field.upper()
dryrun = args.dryrun
prefix = args.prefix
survey = args.survey.lower()
nproc  = args.nproc
realz  = args.realz
'''
config = Configuration(args.config)
config.update_attributes('bound_dist', args)
config.write()
'''
start  = time.time()

# https://www.dur.ac.uk/icc/cosma/cosma5/
fpath  = findfile(ftype='randoms_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix)
opath  = findfile(ftype='randoms_bd', dryrun=dryrun, field=field, survey=survey, prefix=prefix)

if log:
    logfile = findfile(ftype='randoms_bd', dryrun=False, field=field, survey=survey, prefix=prefix, log=True)

    print(f'Logging to {logfile}')

    sys.stdout = open(logfile, 'w')
    
if args.nooverwrite:
    overwrite_check(opath)

call_signature(dryrun, sys.argv)

# Output is sorted by fillfactor.py;   
body      = Table.read(fpath)
boundary  = Table.read(fpath, 'BOUNDARY')

body.sort('CARTESIAN_X')
boundary.sort('CARTESIAN_X')

bids      = boundary['BOUNDID']
boundary  = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]

body      = np.c_[body['CARTESIAN_X'], body['CARTESIAN_Y'], body['CARTESIAN_Z']]

runtime   = calc_runtime(start, 'Reading {:.2f}M randoms'.format(len(body) / 1.e6), xx=body)

split_idx = np.arange(len(body))
split_idx = np.array_split(split_idx, 8 * nproc)

nchunk    = len(split_idx)

runs      = []

for i, idx in enumerate(split_idx):
    split      = body[idx]

    xmin       = split[:,0].min()
    xmax       = split[:,0].max()

    buff       = .1 # Mpc                                                                                                                                                                                 

    # TODO HARDCODE                                                                                                                                                                                        
    complement = (boundary[:,0] > (xmin - 8. - buff)) & (boundary[:,0] < (xmax + 8. + buff))
    complement =  boundary[complement]

    cmin       = complement[:,0].min()
    cmax       = complement[:,0].max()

    print('{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:d}\t{:d}'.format(i, xmin, xmax, cmin, cmax, len(split), len(complement)))

    # leafsize=5                                                                                                                                                                                           
    split      = [x for x in split]
    complement = KDTree(complement)

    runs.append([split, complement])

runtime   = calc_runtime(start, 'Created boundary trees.')

def process_one(run, pid=0):
    split      = run[0]
    complement = run[1]
    '''
    try:
        pid  = multiprocessing.current_process().name.ljust(20)

    except Exception as e:
        print(e)
    '''
    dd, ii     = complement.query(split, k=1)
     
    # del  points

    return  dd.tolist(), ii.tolist()

runtime = calc_runtime(start, 'POOL:  Querying bound dist for body points of {} splits.'.format(nchunk))

now     = time.time()

results = [process_one(runs[0], pid=0)]

split_time  = time.time() - now
split_time /= 60.

runtime = calc_runtime(start, 'POOL:  Expected runtime of {:.3f}.'.format(nchunk * split_time))

with Pool(nproc) as pool:
    for result in tqdm.tqdm(pool.imap(process_one, iterable=runs[1:]), total=len(runs[1:])):
        results.append(result)

    pool.close()
    pool.join()

runtime = calc_runtime(start, 'POOL:  Done with queries')

flat_result = []
flat_ii     = []

for rr in results:
    flat_result   += rr[0]
    flat_ii       += rr[1]

rand               = Table.read(fpath)
rand.sort('CARTESIAN_X')

# print(len(rand))
# print(len(flat_result))

rand['BOUND_DIST'] = np.array(flat_result)
rand['BOUNDID']    = bids[np.array(flat_ii)]

sphere_radius      = rand.meta['RSPHERE']

rand['FILLFACTOR_POISSON'] = rand['FILLFACTOR']
rand['FILLFACTOR'][rand['BOUND_DIST'].data > sphere_radius] = 1.

# CHANGE:  Protect against exactly zero fillfactor (causes division errors). 
rand['FILLFACTOR'] = np.clip(rand['FILLFACTOR'], 1.e-99, None)

runtime = calc_runtime(start, 'Shuffling')

# randomise rows.                                                                                                                                                
idx  = np.arange(len(rand))
idx  = np.random.choice(idx, size=len(idx), replace=False)

rand = rand[idx]

# Bound dist.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query

runtime  = calc_runtime(start, 'Writing {}'.format(opath), xx=rand)

rand.write(opath, format='fits', overwrite=True)

runtime = calc_runtime(start, 'Finished')

if log:
    sys.stdout.close()
