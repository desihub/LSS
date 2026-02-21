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


    

bitmask_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'
bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits'))


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

class get_nobsandmask:
    #doesn't seem like a class should be necessary but couldn't figure out multiprocessing otherwise
    def __init__(self,cat,nproc=128):
        #cat should be an astropy table
        #returns table with same ordering and NOBS_{G,R,Z} and MASKBITS columns
        for col in cat.colnames:
            if not col.isupper(): 
                cat.rename_column(col, col.upper())

        if 'TARGETID' not in cat.colnames:
            cat['TARGETID'] = np.arange(len(cat))

        if 'TARGET_RA' in cat.colnames:
            cat.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])

        if 'INPUT_RA' in cat.colnames:
            cat.rename_columns(['input_ra', 'input_dec'], ['RA', 'DEC'])
        if 'BRICKID' not in cat.colnames:
            from desiutil import brick
            tmp = brick.Bricks(bricksize=0.25)
            cat['BRICKID'] = tmp.brickid(cat['RA'], cat['DEC'])
        # Just some tricks to speed up things up
        self.bid_unique, bidcnts = np.unique(cat['BRICKID'], return_counts=True)
        self.bidcnts = np.insert(bidcnts, 0, 0)
        self.bidcnts = np.cumsum(self.bidcnts)
        self.bidorder = np.argsort(cat['BRICKID'])
        self.cat = cat

        print('adding nobs and mask values to '+str(len(cat))+' rows')
    def wrapper(self,bid_index):

        idx = self.bidorder[self.bidcnts[bid_index]:self.bidcnts[bid_index+1]]
        brickid = self.bid_unique[bid_index]

        ra, dec = self.cat['RA'][idx], self.cat['DEC'][idx]
        tid = self.cat['TARGETID'][idx]
        bitmask,nobsg,nobsr,nobsz = bitmask_radec(brickid, ra, dec)

        data = Table()
        data['idx'] = idx
        data['MASKBITS'] = bitmask
        data['NOBS_G'] = nobsg
        data['NOBS_R'] = nobsr
        data['NOBS_Z'] = nobsz
        data['TARGETID'] = tid

        return data

    # start multiple worker processes
    def get_nobsandmask(self,nproc=128):
        with Pool(processes=nproc) as pool:
            res = pool.map(self.wrapper, np.arange(len(self.bid_unique)))

        res = vstack(res)
        res.sort('idx')
        res.remove_column('idx')
        return res

