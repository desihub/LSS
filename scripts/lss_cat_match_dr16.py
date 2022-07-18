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

parser = argparse.ArgumentParser(description="Match the input catalog to DR16 LSS.")

parser.add_argument("-o", "--out-dir", type=str, required=True, help="Directory to save matched catalog to.")
parser.add_argument("--tracer",help='tracer type to match between',type=str,default='ELG')
parser.add_argument("--version",help='LSS catalog version',type=str,default='test')
parser.add_argument("--specrel",help='LSS catalog version',type=str,default='daily')

args = parser.parse_args()

out_loc = Path(args.out_dir)
if not os.path.isdir(out_loc):
    os.mkdir(out_loc)



if args.specrel == 'daily':
    survey = 'main'
        
if args.specrel == 'guadalupe':
    survey = 'DA02'
if args.specrel  == 'fuji':
    survey = 'SV3'
ROOT = "/global/cfs/cdirs/desi/survey/catalogs/"+survey+"/LSS/"+args.specrel+"/LSScats/"+args.version+"/"
fname = args.tracer+'_full.dat.fits'
with fitsio.FITS(ROOT + fname) as h:
    tab = h[1].read()
    sel = tab['ZWARN'] != 999999 #reject the targets that were not observed
    desi_table = tab[sel]
print("Loaded "+fname+"... matching "+str(len(tab[sel]))+' rows')
    


# Pull out the RA/DEC for use in matching.
desi_ra = desi_table["RA"] 
desi_dec = desi_table["DEC"] 

desi_skycoords = SkyCoord(ra=desi_ra, dec=desi_dec, unit="deg")

# Loads DR16
DR16_ROOT = "/global/cfs/cdirs/sdss/staging/dr16/eboss/lss/catalogs/DR16/"
dr16_fname = "eBOSS_"+args.tracer+"_full_ALLdata-vDR16.fits"

cols_eboss = ["RA", "DEC", "Z", "PLATE", "MJD", "FIBERID","IMATCH"]

with fitsio.FITS(DR16_ROOT + dr16_fname) as h:
    eboss_table = h[1].read_columns(columns=cols_eboss)
    sel = eboss_table['IMATCH'] == 1
    sel |= eboss_table['IMATCH'] == 2
    eboss_table = eboss_table[sel]
print("Loaded "+dr16_fname+"... matching "+str(len(eboss_table))+' rows')
    
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
print("Matched "+str(np.sum(match_keep))+" entries from input catalog to DR16 LSS "+args.tracer+ " catalog.")

# If there are any double matches we'll need to handle that
if np.any(keep_counts) > 1:
    print("Double matches found...")

# Reduces the tables to the matched entries using the indices of matches
desi_keep = Table(desi_table[match_keep])
eboss_keep = Table(eboss_table[idx][match_keep])
eboss_keep.rename_column("Z", "Z_SDSS")
eboss_keep.remove_columns(['RA','DEC'])
joined = hstack([desi_keep, eboss_keep])

# Drops the SDSS RA/DEC from the joined table, since we already have these from
# the DESI portion of the table. 
#del joined["RA"]
#del joined["DEC"]

# Setting the save name.
out_name = args.tracer+"_cat_"+args.specrel+'_'+args.version+"_LSSfull_DR16_match.fits"

joined.write(out_loc / out_name, format="fits", overwrite=True)
    
