# given a full input path, add brick level mask info to it
# Examples:
# 
from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

import LSS.common_tools as common

from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse


time_start = time.time()



parser = argparse.ArgumentParser()
parser.add_argument( '--input_path',required=True)
parser.add_argument('--test', default = 'n', required=False )
parser.add_argument('--RA_col', default = 'RA', required=False )
parser.add_argument('--DEC_col', default = 'DEC', required=False )
parser.add_argument('--n_processes', default = 128, required=False ,type=int)

args = parser.parse_args()

n_processes = args.n_processes


input_path = args.input_path
output_path = args.input_path
    
    

bitmask_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'


def bitmask_radec(brickid, ra, dec):

    brick_index = np.where(bricks['BRICKID']==brickid)[0][0]

    brickname = str(bricks['BRICKNAME'][brick_index])
    if bricks['PHOTSYS'][brick_index]=='N':
        field = 'north'
    elif bricks['PHOTSYS'][brick_index]=='S':
        field = 'south'
    else:
        # raise ValueError
        # Outside DR9 footprint; assign mask bit 16
        bitmask = np.full(len(ra), 2**16, dtype=np.uint8)
        nobsg = 0
        nobsr = 0
        nobsz = 0
        return bitmask,nobsg,nobsr,nobsz

    bitmask_fn = bitmask_dir+'{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)

    bitmask_img = fitsio.read(bitmask_fn)

    header = fits.open(bitmask_fn)[1].header
    w = wcs.WCS(header)

    coadd_x, coadd_y = w.wcs_world2pix(ra, dec, 0)
    coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)

    bitmask = bitmask_img[coadd_y, coadd_x]
    nobsl = []
    bl = ['g','r','z']
    for band in bl:
        nobs_fn = bitmask_dir+'{}/coadd/{}/{}/legacysurvey-{}-nexp-{}.fits.fz'.format(field, brickname[:3], brickname, brickname,band)
        if os.path.isfile(nobs_fn): 
            nobs_img = fitsio.read(nobs_fn)
            nobs = nobs_img[coadd_y,coadd_x]
        else:
            nobs = 0
        nobsl.append(nobs)
    return bitmask,nobsl[0],nobsl[1],nobsl[2]


def wrapper(bid_index):

    idx = bidorder[bidcnts[bid_index]:bidcnts[bid_index+1]]
    brickid = bid_unique[bid_index]

    ra, dec = cat['RA'][idx], cat['DEC'][idx]
    tid = cat['TARGETID'][idx]
    bitmask,nobsg,nobsr,nobsz = bitmask_radec(brickid, ra, dec)

    data = Table()
    data['idx'] = idx
    data['MASKBITS'] = bitmask
    data['NOBS_G'] = nobsg
    data['NOBS_R'] = nobsr
    data['NOBS_Z'] = nobsz
    data['TARGETID'] = tid

    return data


# bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'))
bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits'))


cat = Table(fitsio.read(input_path.replace('global','dvs_ro'))#,columns=[args.RA_col,args.DEC_col]))

if args.RA_col != 'RA':
    cat.rename_columns([args.RA_col],['RA'])
if args.DEC_col != 'DEC':
    cat.rename_columns([args.DEC_col],['DEC'])


print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

if 'TARGETID' not in cat.colnames:
    cat['TARGETID'] = np.arange(len(cat))


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
    res = pool.map(wrapper, np.arange(len(bid_unique)))

res = vstack(res)
res.sort('idx')
res.remove_column('idx')

res.write(output_path,overwrite=True)

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
