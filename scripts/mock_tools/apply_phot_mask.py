from astropy.table import Table
import numpy as np
import time
from multiprocessing import Pool
import fitsio
from astropy import wcs
import os
import gc
from LSS import common_tools as ct
from LSS.globals import main
import argparse
import sys
from astropy.io import fits

ncpus = sys.argv[1]
mocknum = sys.argv[2]
prog = sys.argv[3]

def _run_mask(args):
    fn, nob_fn, idx = args[0], args[1], args[2]
    size = len(idx)
    if os.path.isdir(os.path.dirname(fn)) and all(os.path.exists(nob_fn.format(band)) for band in 'grz'):
    
        bitmask_img, header = fitsio.read(fn, header=True)
        coadd_x, coadd_y = np.round(wcs.WCS(header).wcs_world2pix(cat['RA'][idx], cat['DEC'][idx], 0)).astype(int)
        bitmask = bitmask_img[coadd_y, coadd_x]
        nobs_g, nobs_r, nobs_z = [fitsio.read(nob_fn.format(band))[coadd_y, coadd_x] for band in 'grz']

    else:
        bitmask, nobs_g, nobs_r, nobs_z = [4096]*size, [0]*size, [0]*size, [0]*size
    return bitmask, nobs_g, nobs_r, nobs_z

def get_maskbit_nobs(cat, nproc = int(ncpus), path_to_bricks='/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits', 
                     bitmask_dir = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/'):
    from desiutil import brick
    bricks = Table.read(path_to_bricks)
    #bricks = Catalog.read(path_to_bricks)

    tmp = brick.Bricks(bricksize=0.25)
    cat['BRICKID'] = tmp.brickid(cat['RA'], cat['DEC'])
    cat['BRICKNAME']  = tmp.brickname(cat['RA'], cat['DEC'])
    brickid, bidcnts = np.unique(cat['BRICKID'], return_counts=True)
    bidcnts = np.insert(bidcnts, 0, 0)
    bidcnts = np.cumsum(bidcnts)
    bidorder = np.argsort(cat['BRICKID'])
    bricknames = bricks['BRICKNAME'][np.in1d(bricks['BRICKID'], brickid)]
    regions = bricks['PHOTSYS'][np.in1d(bricks['BRICKID'], brickid)]
    files = [ bitmask_dir+'{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format('south' if reg == 'S' else 'north', brickname[:3], brickname, brickname) for brickname, reg in zip(bricknames, regions) ]
    nobs_fns = [bitmask_dir+'{}/coadd/{}/{}/legacysurvey-{}-nexp-{{}}.fits.fz'.format('south' if reg == 'S' else 'north', brickname[:3], brickname, brickname, None) for brickname, reg in zip(bricknames, regions) ]
    idxs = [bidorder[bidcnts[bid_index]: bidcnts[bid_index+1]] for bid_index in np.arange(len(brickid))]

    print('Initiate Pool with {} processes'.format(nproc), flush=True)
    st=time.time()
    pool = Pool(processes=nproc)
    print('Pool initiated in {} sec'.format(time.time()-st), flush=True)
    
    print('Run maskbit and nobs assignement', flush=True)
    st=time.time()
    res = pool.map(_run_mask, zip(files, nobs_fns, idxs))
    pool.close()
    pool.join()
    cat['MASKBITS'], cat['NOBS_G'], cat['NOBS_R'], cat['NOBS_Z'] = np.hstack([np.vstack(rr) for rr in res])

    gc.collect()
    print('Done in {} sec'.format(time.time()-st), flush=True)
    return cat

if prog=='DARK':
    tag = 'AbacusSummit_v4_1'
    global_type = 'LRG'
elif prog=='BRIGHT':
    tag = 'AbacusSummitBGS_v2'
    global_type = 'BGS'

#file_ = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/forclustering/cutsky_abacusHF_DR2_ELG_z0p950_zcut_0p8to1p6_clustering.dat.fits'
##file_ = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_nomask_%d_v2.fits' % int(mocknum)
##TMEPfile_ = '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/{TAG}/forFA{MOCK}_nomask.fits'.format(MOCK=mocknum, TAG=tag) 
file_ = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/forFA202_noimagingmask.fits'
#file_ = '/dvs_ro/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/testcontaminants.fits'

#file_ = '/pscratch/sd/e/efdez/Uchuu-GLAM/GLAM/mocks_altmtl/LRG/GLAM-Uchuu_LRG_100_Y3_cut_sky_clustering.h5'
#file_ = '/pscratch/sd/e/efdez/Uchuu/LSS/scripts/mock_tools/DA2/LRG_NGC_12_clustering.ran.fits'
#/pscratch/sd/z/zxzhai/DESI/PreMocks/SecondGenMocks/AbacusSummit_v4_1/forFA'+str(mocknum)+'_nomasking.fits'

'''
import h5py
import hdf5plugin #need to be in the cosmodesi test environment, as of Sep 4th 25
cat = Table()
with h5py.File(file_) as fn:
    columns = fn.keys()
    for col in columns:
        cat[col] = fn[col][:]
'''
cat = Table.read(file_)



print('size before cutting photmask is', len(cat))
cat = get_maskbit_nobs(cat)

print('size after matching photmask is', len(cat))

mainp = main(tp = global_type, specver = 'loa-v1')
cat = ct.cutphotmask(cat, bits=mainp.imbits)

st=time.time()
print('size after cutting photmask is', len(cat))
#TEMPout_file_name='/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/{TAG}/forFA{MOCK}.fits'.format(MOCK=mocknum, TAG=tag)
#out_file_name='/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_%d_v2.fits' % int(mocknum)
#out_file_name = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/forclustering/masked_cutsky_abacusHF_DR2_ELG_z0p950_zcut_0p8to1p6_clustering.dat.fits'
##out_file_name = 'test_glam_lrg.fits'
out_file_name = '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/forFA202.fits'
#out_file_name = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/forFA_QSO.fits'
ct.write_LSS_scratchcp(cat, out_file_name, extname='TARGETS')
print('Done writing in {} sec'.format(time.time()-st), flush=True)
fits.setval(out_file_name, 'OBSCON', value=prog, ext=1)

