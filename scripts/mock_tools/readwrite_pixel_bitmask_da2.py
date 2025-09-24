# Get bitmask values from pixel-level per-brick masks for a catalog
# Examples:
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer lrg --input catalog.fits --output catalog_lrgmask_v1.1.npy
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer lrg --input /global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits --output $CSCRATCH/temp/randoms-1-0-lrgmask_v1.1.fits
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer lrg --input /global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_south_1.0.0_basic.fits --output $CSCRATCH/temp/dr9_lrg_south_1.0.0_lrgmask_v1.1.fits

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


time_start = time.time()


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tracer', required=True)
parser.add_argument('-i', '--input', required=True)
parser.add_argument('-o', '--output', required=False)
parser.add_argument('-v', '--version', default='none', required=False)
parser.add_argument('-rv', '--tarver', default='targetsDR9v1.1.1', required=False)
parser.add_argument( '--cat_type', default='targets')#, choices=['targets','ran','obielg','Ab2ndgen'],required=False)
parser.add_argument( '--reg', default='north', choices=['north','south'],required=False)
parser.add_argument( '--secgen_ver', default=None, required=False)
parser.add_argument( '--nproc', default=128, required=False)

args = parser.parse_args()

n_processes = int(args.nproc)

input_path = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+args.input+args.tarver+'.fits'
output_path = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+args.input+args.tarver+'_'+args.tracer+'imask.fits'
if args.cat_type == 'ran':
    input_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-'+str(args.input)+'.fits'
    output_path = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-'+str(args.input)+args.tracer+'imask.fits'
if args.cat_type == 'obielg':
    input_path = '/global/cfs/cdirs/desi/survey/catalogs/image_simulations/ELG/dr9/Y1/'+args.reg+'/file0_rs0_skip0/merged/matched_input_full.fits'
    output_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/elg_obiwan_'+args.reg+'_matched_input_full_'+args.tracer+'_imask.fits'
if args.cat_type == 'Ab2ndgen':
    if args.secgen_ver is None:
        mockdir = '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit/'
    else:
        mockdir = '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/%s/' %  args.secgen_ver
    input_path = mockdir+'forFA'+args.input+'.fits'
    output_path = mockdir+'forFA'+args.input+'_matched_input_full_'+args.tracer+'_imask.fits'
if args.cat_type == 'Y1EZmock':
    mockdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/forFA/'
    input_path = mockdir+'forFA'+args.input+'.fits'
    output_path = mockdir+'forFA'+args.input+'_matched_input_full_'+args.tracer+'_imask.fits'

if args.cat_type == 'Generic':
    input_path = args.input
    output_path = args.output 


tracer = args.tracer.lower()
version = args.version

version_dict = {'lrg': 'v1.1', 'elg': 'v1'}
if version=='none':
    version = version_dict[tracer]

bitmask_dir = '/global/cfs/cdirs/desi/survey/catalogs/brickmasks/{}/{}'.format(tracer.upper(), version)

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
bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits'))

#try:
#    cat = Table(fitsio.read(input_path, rows=None, columns=['RA', 'DEC', 'BRICKID', 'TARGETID']))
#except ValueError:
#    cat = Table(fitsio.read(input_path, rows=None, columns=['RA', 'DEC', 'TARGETID']))

if args.cat_type == 'obielg':
    cat = Table(fitsio.read(input_path,columns=['input_ra','input_dec']))
    cat.rename_columns(['input_ra', 'input_dec'], ['RA', 'DEC'])
else:
    cat = Table(fitsio.read(input_path))

print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

if 'TARGETID' not in cat.colnames:
    cat['TARGETID'] = np.arange(len(cat))

if 'TARGET_RA' in cat.colnames:
    cat.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])

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
    res.write(output_path)
else:
    np.write(output_path, np.array(res['{}_mask'.format(tracer)]))

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
