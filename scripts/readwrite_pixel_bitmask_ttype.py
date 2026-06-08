# Get bitmask values from pixel-level per-brick masks defined for a particular tracer type (lrg and elg available)
# Examples:
# srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python readwrite_pixel_bitmask_ttype.py --tracer lrg --input_fn catalog.fits 

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse

import LSS.common_tools as common


time_start = time.time()

n_processes = 128

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tracer',help='the tracer type for the mask', required=True)
parser.add_argument('-i', '--input_fn',help='the input file name, without .fits', required=True)
parser.add_argument('-id', '--input_dir',help='the input directory',default='')
parser.add_argument('--output_fn',help='the output file name, by default will match input with {tracer}imask.fits at the end', default=None)
parser.add_argument('--output_dir',help='the output directory',default='')
parser.add_argument('-v', '--version', default='none', required=False)
parser.add_argument('--ra_col',default='RA')
parser.add_argument('--dec_col',default='DEC')
args = parser.parse_args()


input_path = args.input_dir+'/'+args.input_fn+'.fits'
if args.output_fn is None:
    out_fn = args.input_fn +'_'+args.tracer+'imask.fits'
else:
    out_fn = args.output_fn
output_path = args.output_dir+'/'+out_fn
   


tracer = args.tracer.lower()
version = args.version

version_dict = {'lrg': 'v1.1', 'elg': 'v1'}
if version=='none':
    version = version_dict[tracer]

bitmask_dir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/brickmasks/{}/{}'.format(tracer.upper(), version)

# input_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits'
# output_path = '/global/cscratch1/sd/rongpu/temp/randoms-1-0-lrgmask_v1.fits'

if os.path.isfile(output_path):
    raise ValueError(output_path+' already exists!')


def bitmask_radec(brickid, ra, dec):

    brick_index = np.where(bricks['BRICKID']==brickid)[0][0]

    brickname = str(bricks['BRICKNAME'][brick_index])
    if bricks['PHOTSYS'][brick_index]=='N':
        field = 'north'
    elif bricks['PHOTSYS'][brick_index]=='S':
        field = 'south'
    else:
        # raise ValueError
        # Outside DR9 footprint; assign mask bit 7
        bitmask = np.full(len(ra), 2**7, dtype=np.uint8)
        return bitmask

    # bitmask_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)
    bitmask_fn = os.path.join(bitmask_dir, '{}/coadd/{}/{}/{}-{}mask.fits.gz'.format(field, brickname[:3], brickname, brickname, tracer))

    bitmask_img = fitsio.read(bitmask_fn)

    header = fits.open(bitmask_fn)[1].header
    w = wcs.WCS(header)

    coadd_x, coadd_y = w.wcs_world2pix(ra, dec, 0)
    coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)

    bitmask = bitmask_img[coadd_y, coadd_x]

    return bitmask


def wrapper(bid_index):

    idx = bidorder[bidcnts[bid_index]:bidcnts[bid_index+1]]
    brickid = bid_unique[bid_index]

    ra, dec = cat['RA'][idx], cat['DEC'][idx]
    tid = cat['TARGETID'][idx]
    bitmask = bitmask_radec(brickid, ra, dec)

    data = Table()
    data['idx'] = idx
    data['{}_mask'.format(tracer)] = bitmask
    data['TARGETID'] = tid

    return data


# bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'))
bricks = Table(fitsio.read('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits'))


cat1row = fitsio.read(input_path,rows=1)
cols = [args.ra_col,args.dec_col,'TARGETID','BRICKID']
cols2read = []
cols_in_f = list(cat1row.dtype.names)
for col in cols:
    if col in cols_in_f:
        cols2read.append(col)
    if col.lower() in cols_in_f:
        cols2read.append(col)

del cat1row

cat = Table(fitsio.read(input_path.replace('global','dvs_ro'),columns=cols2read))

print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

if 'TARGETID' not in cat.colnames:
    cat['TARGETID'] = np.arange(len(cat))


cat.rename_columns([args.ra_col, args.dec_col], ['RA', 'DEC'])

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

if output_path.endswith('.fits'):
    common.write_LSS_scratchcp(res,output_path)
else:
    np.write(output_path, np.array(res['{}_mask'.format(tracer)]))

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
