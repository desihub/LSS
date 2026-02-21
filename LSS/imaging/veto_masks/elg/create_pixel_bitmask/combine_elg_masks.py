# Create the combined ELG mask from the individual masks

# To run on a single interactive node:
# python combine_elg_masks.py "south 1 0"

# Full production run on interactive nodes:
# salloc -N 3 -C haswell -q interactive -t 04:00:00
# srun --wait=0 --ntasks-per-node 1 payload.sh tasks.txt ; exit

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool
import argparse


ts_bit = 0
custom_bit = 1
gaia_bit = 2
elg_custom_bit = 3

ts_maskbits = [1, 12, 13]  # DESI targeting mask bits

custom_input_dir = '/global/cfs/cdirs/desi/users/rongpu/data/veto_masks/all_tracers_custom_mask/v1'
elg_custom_input_dir = '/global/cscratch1/sd/rongpu/desi/veto_masks/elg/dev/elg_custom/v1'
gaia_input_dir = '/global/cscratch1/sd/rongpu/desi/veto_masks/elg/dev/gaiamask/v1'
output_dir = '/global/cscratch1/sd/rongpu/desi/veto_masks/elg/v1'
# output_dir = '/global/cfs/cdirs/desi/survey/catalogs/brickmasks/ELG/v1'

parser = argparse.ArgumentParser()
parser.add_argument('args')
args = parser.parse_args()
field, n_task, task_id = args.args.split()
n_task, task_id = int(n_task), int(task_id)

n_processes = 32


################################
debug = False
################################


def get_combined_bitmask(brick_index):

    brickname = str(bricks['brickname'][brick_index])

    output_path = os.path.join(output_dir, '{}/coadd/{}/{}/{}-elgmask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    if os.path.isfile(output_path):
        return None

    # print(output_path)

    img_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)

    header = fitsio.read_header(img_fn, ext=1)
    header_keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    hdict = {}
    for keyword in header_keywords:
        hdict[keyword] = header[keyword]
    if (hdict['CTYPE1']!='RA---TAN') or (hdict['CTYPE2']!='DEC--TAN'):
        raise ValueError

    bitmask = np.full((3600, 3600), 0, dtype=np.uint8)

    # targeting maskbits
    maskbits = fitsio.read(img_fn, ext='MASKBITS')
    mask_ts = np.full(maskbits.shape, False)
    for bit in ts_maskbits:
        mask_ts |= (maskbits & 2**bit)>0
    bitmask[mask_ts] += 2**ts_bit

    customm_path = os.path.join(custom_input_dir, 'coadd/{}/{}/{}-custommask.fits.gz'.format(brickname[:3], brickname, brickname))
    if os.path.isfile(customm_path):  # there is no custom masking if the file does not exist
        customm = fitsio.read(customm_path).astype(np.uint8)
        bitmask[customm!=0] += 2**custom_bit

    elg_customm_path = os.path.join(elg_custom_input_dir, 'coadd/{}/{}/{}-elg-custommask.fits.gz'.format(brickname[:3], brickname, brickname))
    if os.path.isfile(elg_customm_path):  # there is no elg_custom masking if the file does not exist
        elg_customm = fitsio.read(elg_customm_path).astype(np.uint8)
        bitmask[elg_customm!=0] += 2**elg_custom_bit

    gaiam_path = os.path.join(gaia_input_dir, '{}/coadd/{}/{}/{}-gaiamask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    gaiam = fitsio.read(gaiam_path).astype(np.uint8)
    bitmask[gaiam!=0] += 2**gaia_bit

    if not os.path.isdir(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except:
            pass

    fitsio.write(output_path, bitmask, compress='GZIP', header=hdict, clobber=True)

    # return bitmask
    return None


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
    res = pool.map(get_combined_bitmask, bricks_list, chunksize=1)

print('combine_elg_mask {} {} {} Done!'.format(field, n_task, task_id), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

