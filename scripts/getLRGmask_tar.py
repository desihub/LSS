# Get LRG bitmasks for a catalog
# originally written by Rongpu Zhou
# Examples:
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --input catalog.fits --output catalog_lrgmask.npy
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --input /global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits --output $CSCRATCH/temp/randoms-1-0-lrgmask_v1.fits
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --input /global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_south_1.0.0_basic.fits --output $CSCRATCH/temp/dr9_lrg_south_1.0.0_lrgmask_v1.fits

from __future__ import division, print_function
from functools import partial
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse


time_start = time.time()

#bitmask_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/v1'
bitmask_dir = '/global/cfs/cdirs/desi/survey/catalogs/brickmasks/LRG/v1'

n_processes = 32

##################
debug = False
##################

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--survey", help="e.g., SV3 or main",default='main')

#parser.add_argument('-i', '--input', required=True)
#parser.add_argument('-o', '--output', required=True)
args = parser.parse_args()

lssdir = args.basedir +'/'+args.survey+'/LSS/'

tp = 'LRG'



def bitmask_radec(brickid, ra, dec):

    brick_index = np.where(bricks['BRICKID']==brickid)[0][0]

    brickname = str(bricks['BRICKNAME'][brick_index])
    if bricks['PHOTSYS'][brick_index]=='N':
        field = 'north'
    elif bricks['PHOTSYS'][brick_index]=='S':
        field = 'south'
    else:
        raise ValueError
    # bitmask_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)
    bitmask_fn = os.path.join(bitmask_dir, '{}/coadd/{}/{}/{}-lrgmask.fits.gz'.format(field, brickname[:3], brickname, brickname))

    bitmask_img = fitsio.read(bitmask_fn)

    header = fits.open(bitmask_fn)[1].header
    w = wcs.WCS(header)

    coadd_x, coadd_y = w.wcs_world2pix(ra, dec, 0)
    coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)

    bitmask = bitmask_img[coadd_y, coadd_x]

    return bitmask

def wrapper(bid_index,bidorder,bidcnts,bid_unique,cat):

    idx = bidorder[bidcnts[bid_index]:bidcnts[bid_index+1]]
    brickid = bid_unique[bid_index]

    ra, dec = cat['RA'][idx], cat['DEC'][idx]

    bitmask = bitmask_radec(brickid, ra, dec)

    data = Table()
    data['idx'] = idx
    data['lrg_mask'] = bitmask
    data['TARGETID'] = cat['TARGETID'][idx]

    return data


def mkfile(input_path,output_path):
    try:
        cat = fitsio.read(input_path, rows=None, columns=['lrg_mask'])
        return 'file already has lrg_mask column'
    except:
        print('adding lrg_mask column') 

    try:
        cat = Table(fitsio.read(input_path, rows=None, columns=['RA', 'DEC', 'BRICKID','TARGETID']))
    except ValueError:
        cat = Table(fitsio.read(input_path, rows=None, columns=['RA', 'DEC','TARGETID']))

    print(len(cat))

    #for col in cat.colnames:
    #    cat.rename_column(col, col.upper())

    #if 'TARGET_RA' in cat.colnames:
    #    cat.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])

    if 'BRICKID' not in cat.colnames:
        from desiutil import brick
        tmp = brick.Bricks(bricksize=0.25)
        cat['BRICKID'] = tmp.brickid(cat['RA'], cat['DEC'])

    # Just some tricks to speed up things up
    bid_unique, bidcnts = np.unique(cat['BRICKID'], return_counts=True)
    bidcnts = np.insert(bidcnts, 0, 0)
    bidcnts = np.cumsum(bidcnts)
    bidorder = np.argsort(cat['BRICKID'])


    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(partial(wrapper,bidorder=bidorder,bidcnts=bidcnts,bid_unique=bid_unique,cat=cat), np.arange(len(bid_unique)))
        #partial(func, b=second_arg), a_args

    res = vstack(res)
    res.sort('idx')
    res.remove_column('idx')

    cat = Table(fitsio.read(input_path))

    if len(cat) != len(res):
        print('mismatched lengths, somehow get brick mask removed data!!!')

    else:
        res = join(cat,res,keys=['TARGETID'])
        if output_path.endswith('.fits'):
            res.write(output_path,overwrite=True)
        else:
            np.write(output_path, np.array(res['lrg_mask']))
        del cat
        del res
        print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

# bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'))
bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits'))

if debug:
    rows = np.arange(int(1e3))
else:
    rows = None

if args.survey == 'main':
    input_path = lssdir+tp+'targetsDR9v1.1.1.fits'
    output_path = input_path #we will over-write, just adding new column

mkfile(input_path,output_path)


