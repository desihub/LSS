# Create pixel-level per-brick LRG GAIA mask

# To run on a single interactive node (for debugging):
# python lrg_gaiamask.py "south 1 0"

# Full production run on interactive nodes:
# salloc -N 3 -C haswell -q interactive -t 04:00:00
# srun --wait=0 --ntasks-per-node 1 payload.sh tasks.txt ; exit

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from astropy.io import fits
from astropy import wcs

from scipy.interpolate import RectBivariateSpline

from astropy import units as u
from astropy.coordinates import SkyCoord

from multiprocessing import Pool
import argparse


gaia_bit = 0

gaia_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_lrg_mask_v1.fits'

output_dir = '/global/cscratch1/sd/rongpu/desi/veto_masks/lrg/dev/gaiamask/v1'

parser = argparse.ArgumentParser()
parser.add_argument('args')
args = parser.parse_args()
field, n_task, task_id = args.args.split()
n_task, task_id = int(n_task), int(task_id)

n_processes = 32

gaia_columns = ['RA', 'DEC', 'mask_mag', 'radius_'+field]

################################
debug = False
################################


def get_pixel_bitmask(brick_index):

    brickname = str(bricks['brickname'][brick_index])

    output_path = os.path.join(output_dir, '{}/coadd/{}/{}/{}-gaiamask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    if os.path.isfile(output_path):
        return None

    # print(output_path)

    ra1, dec1 = [bricks['ra'][brick_index]], [bricks['dec'][brick_index]]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')

    search_radius = stars['radius'].max() + 0.2*3600
    _, idx2, d2d, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    # print(len(idx2))

    d2d = np.array(d2d.to(u.arcsec))
    mask = d2d < (stars['radius'][idx2] + 0.2*3600)
    # print(np.sum(mask))
    idx2 = idx2[mask]
    d2d = d2d[mask]

    sky2_brick = SkyCoord(ra2[idx2]*u.degree, dec2[idx2]*u.degree, frame='icrs')
    stars_brick = stars[idx2].copy()

    img_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)

    header = fits.open(img_fn)[1].header
    header_keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    hdict = {}
    for keyword in header_keywords:
        hdict[keyword] = header[keyword]
    if (hdict['CTYPE1']!='RA---TAN') or (hdict['CTYPE2']!='DEC--TAN'):
        raise ValueError

    w = wcs.WCS(header)
    naxis1 = header['NAXIS1']  # Length of the *second* index of the 2-D array
    naxis2 = header['NAXIS2']  # Length of the *first* index of the 2-D array

    if naxis1!=3600 or naxis2!=3600:
        raise ValueError

    binsize = 1000
    pix_x_spline, pix_y_spline = np.arange(-binsize, naxis1+2*binsize, binsize), np.arange(-binsize, naxis2+2*binsize, binsize)
    xx, yy = np.meshgrid(pix_x_spline, pix_y_spline)
    pix_ra_spline, pix_dec_spline = w.wcs_pix2world(xx, yy, 0)

    # Get rid of the discontinuity in RA for bricks at touching RA=0
    if bricks['ra1'][brick_index]==0:
        # wraparound to RA=0
        pix_ra_spline = ((pix_ra_spline+180.)%360.)-180.
    elif bricks['ra2'][brick_index]==360:
        # wraparound to RA=360
        pix_ra_spline = (pix_ra_spline-180.)%360.+180.

    interp_ra = RectBivariateSpline(pix_y_spline, pix_x_spline, pix_ra_spline)
    interp_dec = RectBivariateSpline(pix_y_spline, pix_x_spline, pix_dec_spline)

    chunk_size = 400
    if naxis1%chunk_size!=0:
        raise ValueError
    n_chunks = naxis1//chunk_size

    bitmask_i = []
    for i in range(n_chunks):

        bitmask_j = []
        for j in range(n_chunks):

            bitmask = np.full((chunk_size, chunk_size), 0, dtype=np.uint8)

            pix_ra = interp_ra(i*chunk_size+np.arange(chunk_size), j*chunk_size+np.arange(chunk_size)).flatten()
            pix_dec = interp_dec(i*chunk_size+np.arange(chunk_size), j*chunk_size+np.arange(chunk_size)).flatten()

            sky1_brick = SkyCoord([np.mean(pix_ra)]*u.degree, [np.mean(pix_dec)]*u.degree, frame='icrs')

            _, idx2, d2d, _ = sky2_brick.search_around_sky(sky1_brick, seplimit=search_radius*u.arcsec)

            d2d = np.array(d2d.to(u.arcsec))
            mask = d2d < (stars_brick['radius'][idx2] + 0.2*chunk_size/3600*3600)
            idx2 = idx2[mask]
            # print(len(idx2))

            if len(idx2)==0:
                bitmask_j.append(bitmask)
                continue

            # RA/DEC to unit cartesian vectors
            pix_cx = np.cos(np.radians(pix_ra))*np.cos(np.radians(pix_dec))
            pix_cy = np.sin(np.radians(pix_ra))*np.cos(np.radians(pix_dec))
            pix_cz = np.sin(np.radians(pix_dec))

            star_ra = stars_brick['RA'][idx2]
            star_dec = stars_brick['DEC'][idx2]
            mask_radii = stars_brick['radius'][idx2]

            # RA/DEC to unit cartesian vectors
            star_cx = np.cos(np.radians(star_ra))*np.cos(np.radians(star_dec))
            star_cy = np.sin(np.radians(star_ra))*np.cos(np.radians(star_dec))
            star_cz = np.sin(np.radians(star_dec))

            mat1 = np.array([pix_cx, pix_cy, pix_cz]).T
            mat2 = np.array([star_cx, star_cy, star_cz])
            del pix_cx, pix_cy, pix_cz, star_cx, star_cy, star_cz

            dist = np.dot(mat1, mat2)
            mask = np.any(dist>np.cos(np.radians(mask_radii/3600.)), axis=1).reshape(chunk_size, chunk_size)
            bitmask[mask] += 2**gaia_bit

            bitmask_j.append(bitmask)

        bitmask_i.append(np.hstack(bitmask_j))

    bitmask = np.vstack(bitmask_i)

    if not os.path.isdir(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except:
            pass
    fitsio.write(output_path, bitmask, compress='GZIP', header=hdict, clobber=True)

    # return bitmask
    return None


hdu = fits.open(gaia_path)
stars = Table()
for col in gaia_columns:
    stars[col] = np.copy(hdu[1].data[col])
# print(len(stars))

stars.rename_column('radius_'+field, 'radius')

##############################################################

ra2, dec2 = np.array(stars['RA']), np.array(stars['DEC'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

# Create KD tree cache
ra0, dec0 = [0.], [0.]
sky0 = SkyCoord(ra0*u.degree, dec0*u.degree, frame='icrs')
search_radius = stars['radius'].max() + 0.2*3600
_, _, _, _ = sky2.search_around_sky(sky0, seplimit=search_radius*u.arcsec)

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/survey-bricks-dr9-{}.fits.gz'.format(field, field)))
# print(len(bricks))

# ########################## single brick ##########################
# mask = bricks['brickname']=='1092p320'
# brick_index = np.where(mask)[0][0]
# bitmask = get_pixel_bitmask(brick_index)
# ##################################################################

###################################################
if debug:
    np.random.seed(213)
    idx = np.random.choice(len(bricks), size=(10000), replace=False)
    bricks = bricks[idx]
###################################################

# random shuffle
np.random.seed(4891)
bricks_list = np.random.choice(len(bricks), size=len(bricks), replace=False)
# split among the Cori nodes
bricks_list_split = np.array_split(bricks_list, n_task)
bricks_list = bricks_list_split[task_id]
print('Number of bricks in this node:', len(bricks_list))

time_start = time.time()
with Pool(processes=n_processes) as pool:
    res = pool.map(get_pixel_bitmask, bricks_list, chunksize=1)

print('lrg_gaiamask {} {} {} Done!'.format(field, n_task, task_id), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

