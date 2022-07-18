#!/usr/bin/env python3

import argparse
import glob
import os
from pathlib import Path

from astropy import units as u
from astropy.coordinates import match_coordinates_sky, SkyCoord
from astropy.table import Table, hstack
import numpy as np
import fitsio

parser = argparse.ArgumentParser(description="Match the input catalog to DR16Q.")

parser.add_argument("-o", "--out-dir", type=str, required=True, help="Directory to save matched catalog to.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-f", "--fuji", action="store_true", help="Match against Fuji catalog.")
group.add_argument("-g", "--guadalupe", action="store_true", help="Match against Guadalupe catalog.")
group.add_argument("-d", "--daily", action="store_true", help="Match against daily catalog.")
group.add_argument("-a", "--all", action="store_true", help="Match against combined catalogs.")

args = parser.parse_args()

out_loc = Path(args.out_dir)
if not os.path.isdir(out_loc):
    os.mkdir(out_loc)

# Load each of the two releases individually
releases = []
if args.guadalupe or args.all:
    releases.append("guadalupe")
if args.fuji or args.all:
    releases.append("fuji")
if args.daily or args.all:
    releases.append("daily")

desi_tables = {}

for r in releases:
    if r == 'daily':
        ROOT = "/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/"
        fname = "QSO_catalog.fits"
    else:
        ROOT = f"/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/{r}/"
        fname = f"QSO_cat_{r}_healpix.fits"

    with fitsio.FITS(ROOT + fname) as h:
        desi_tables[r] = h[1].read()
    print(f"Loaded {fname}...")

# Combine the two releases into a single table
desi_table_combined = desi_table = np.concatenate([desi_tables[k] for k in releases])

# Pull out the RA/DEC for use in matching.
desi_ra = np.asarray([i["TARGET_RA"] for i in desi_table])
desi_dec = np.asarray([i["TARGET_DEC"] for i in desi_table])

desi_skycoords = SkyCoord(ra=desi_ra, dec=desi_dec, unit="deg")

# Loads DR16
DR16Q_ROOT = "/global/cfs/cdirs/sdss/staging/dr16/eboss/qso/DR16Q/"
dr16q_fname = "DR16Q_v4.fits"

cols_eboss = ["RA", "DEC", "Z", "PLATE", "MJD", "FIBERID"]

with fitsio.FITS(DR16Q_ROOT + dr16q_fname) as h:
    eboss_table = h[1].read_columns(columns=cols_eboss)
    print(f"Loaded {dr16q_fname}...")
    
eboss_ra = np.asarray([i["RA"] for i in eboss_table])
eboss_dec = np.asarray([i["DEC"] for i in eboss_table])
eboss_skycoords = SkyCoord(ra=eboss_ra, dec=eboss_dec, unit="deg")

# This is the line that actually matches the two table RA/DECs to each other
print("Matching...")
idx, sep2d, dist3d = match_coordinates_sky(desi_skycoords, eboss_skycoords)

# 2d seperation in arc seconds to constrain our search radius.
d2d = np.asarray(sep2d.to(u.arcsec))

# Keep everything whose match is within 1 arcsecond
# Eseentially deciding everything that close is "correct"
match_keep = d2d < 1
_, keep_counts = np.unique(idx[match_keep], return_counts=True)
print(f"Matched {np.sum(match_keep)} entries from input catalog to DR16Q.")

# If there are any double matches we'll need to handle that
if np.any(keep_counts) > 1:
    print("Double matches found...")

# Reduces the tables to the matched entries using the indices of matches
desi_keep = Table(desi_table[match_keep])
eboss_keep = Table(eboss_table[idx][match_keep])
eboss_keep.rename_column("Z", "Z_SDSS")
joined = hstack([desi_keep, eboss_keep])

# Drops the SDSS RA/DEC from the joined table, since we already have these from
# the DESI portion of the table. 
del joined["RA"]
del joined["DEC"]

# Setting the save name.
out_name = "QSO_cat_fujilupe_healpix_DR16Q_match.fits"
if args.fuji:
    out_name = "QSO_cat_fuji_healpix_DR16Q_match.fits"
elif args.guadalupe:
    out_name = "QSO_cat_guadalupe_healpix_DR16Q_match.fits"
elif args.daily:
    out_name = "QSO_cat_daily_tile_DR16Q_match.fits"

joined.write(out_loc / out_name, format="fits", overwrite=True)
    
