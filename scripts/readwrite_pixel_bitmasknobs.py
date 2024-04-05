# Get bitmask values from pixel-level per-brick masks for a catalog
# Examples:
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python read_pixel_bitmaskmasknobs.py 
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
parser.add_argument( '--cat_type', default='obielg', choices=['obielg', 'abacus'],required=False)
parser.add_argument( '--reg', default='north', choices=['north','south'],required=False)
parser.add_argument('--abacus_fa_num', default = 0, required = False)
parser.add_argument('--do_randoms', default = 'n', choices = ['n','y'], required = False)
parser.add_argument('--random_tracer', default = 'LRG', required = False)
parser.add_argument('--mock_number', default = 0, required = False)
parser.add_argument('--outdir', default = '', required=False )
parser.add_argument('--overwrite', default = 'n', required=False )
parser.add_argument('--test', default = 'n', required=False )
parser.add_argument('--n_processes', default = 128, required=False ,type=int)

args = parser.parse_args()

n_processes = args.n_processes

if args.cat_type == 'obielg':
    input_path = '/global/cfs/cdirs/desi/survey/catalogs/image_simulations/ELG/dr9/Y1/'+args.reg+'/file0_rs0_skip0/merged/matched_input_full.fits'
    output_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/elg_obiwan_'+args.reg+'_matched_input_full_masknobs.fits'
    
if args.cat_type == 'abacus':
    if args.do_randoms == 'n':
        fa_num = args.abacus_fa_num
        str_fa_num = str(fa_num)
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/"
        input_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/forFA' + str_fa_num + '.fits'
        if args.outdir != '':
            output_path = args.outdir + "/" + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"
        elif args.outdir == '':
            output_path = input_dir + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"
            print(output_path)
            
    elif args.do_randoms == 'y':
        ran_tr = args.random_tracer
        mockno = args.mock_number
        print("Running for Mock %s on Tracer %s"%(mockno, ran_tr))
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/"%(mockno)
        input_path = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/%s_1_full_noveto.ran.fits"%(mockno, ran_tr)
        if args.outdir != '':
            output_path = args.outdir + "/" + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"
            print("Output to " + output_path)
        elif args.outdir == '':
            output_path = input_dir + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"
            print("Output to " + output_path)
    

bitmask_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'

# input_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits'
# output_path = '/global/cscratch1/sd/rongpu/temp/randoms-1-0-lrgmask_v1.fits'
fe = False
if os.path.isfile(output_path):
    fe = True
    if args.overwrite == 'n' and args.test == 'n':
        raise ValueError(output_path+' already exists!')
    if args.overwrite == 'n' and args.test == 'y':
        print('will run and get timing, but no output will be written (because path exists)')
    if args.overwrite == 'y':
        print('will overwrite '+output_path)


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

#if output_path.endswith('.fits'):
#ow = False
if fe == False or args.overwrite == 'y':
    #ow = True
    res.write(output_path,overwrite=True)
#else:
#    np.write(output_path, np.array(res['masknobs']))

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
