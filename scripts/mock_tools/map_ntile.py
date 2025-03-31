#In order to run this script, you need the cosmodesi environment
#source /global/common/software/desi/desi_environment.sh main
#source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#Then, request two nodes interactively (you need 2 for memory, the code works with ~400GB of data)
#salloc -N 1 --qos interactive --time 04:00:00 --constraint cpu --account desi 
#then, run
#to run this, just python /global/cfs/cdirs/desi/users/schiarenza/LSS/scripts/mock_tools/map_ntile.py
import os
import argparse
import numpy as np
import healpy as hp
import fitsio
import tqdm

def radec2thphi(ra, dec):
    return (-dec + 90.) * np.pi / 180., ra * np.pi / 180.

def map_ntile(ra, dec, ntile, nside=1024):
    print("Starting map_ntile function...")  # Debugging start message

    # Convert RA, DEC to theta, phi
    theta, phi = radec2thphi(ra, dec)
    # Map to healpix pixels
    pix = hp.ang2pix(nside, theta, phi, nest=True)

    map_c = np.zeros(12 * nside**2)
    map_ntile = np.zeros(12 * nside**2)

    print(f"Processing {len(pix)} objects...")  # Show total objects to process

    # Use tqdm for progress tracking
    for ii in tqdm.tqdm(range(len(pix)), desc="Mapping NTILE"):
        pi = pix[ii]
        nt = ntile[ii]
        map_c[pi] += 1.
        map_ntile[pi] += nt

    # Normalize NTILE map
    sel = map_c > 0
    map_ntile[sel] /= map_c[sel]

    print(f"NTILE map computed. Min: {np.min(map_ntile)}, Max: {np.max(map_ntile)}")

    return map_ntile

# Argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--prog', help='tracer to be selected', type=str, choices=['dark', 'bright'], default='dark')
parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main', 'Y1'], default='Y1')
parser.add_argument('--verspec', help='version for redshifts', type=str, default='iron')
parser.add_argument('--version', help='catalog version', type=str, default='v1.5')
parser.add_argument('--nside', help='nside for healpix map', type=int, default=1024)
parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=18)
args = parser.parse_args()

# Set input directory
indir = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.verspec}/LSScats/{args.version}/'

# Determine tracer
tracer = 'QSO' if args.prog == 'dark' else 'BGS_BRIGHT'

nside = args.nside
print(f"Using nside {nside}")

# Read the random files
if args.nran == 1:
    ran_fn = f'{indir}/{tracer}_1_full.ran.fits'
    ran_cat = fitsio.read(ran_fn)
else:
    ran_fns = [f'{indir}/{tracer}_{idx}_full.ran.fits' for idx in range(args.nran)]
    print("Loading randoms...")

    ra_list = []
    dec_list = []
    ntile_list = []

    for fn in ran_fns:
        print(f"Reading {fn}...")
        # Read the columns directly as numpy arrays
        data = fitsio.read(fn)
        ra_list.append(data['RA'])
        dec_list.append(data['DEC'])
        ntile_list.append(data['NTILE'])

    # Concatenate the columns
    ra_concat = np.concatenate(ra_list)
    dec_concat = np.concatenate(dec_list)
    ntile_concat = np.concatenate(ntile_list)

    print(f"All {args.nran} randoms loaded!")

# Run the map_ntile function
map_nt = map_ntile(ra_concat, dec_concat, ntile_concat, nside=nside)

# Save the output
save_fn = f'{indir}/healpix_map_ntile_{args.prog}_nside_{args.nside}.fits'
hp.write_map(save_fn, map_nt, overwrite=True)
print(f"NTILE map saved to {save_fn}")
