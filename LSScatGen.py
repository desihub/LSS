# This script will swipe through the output of the data reduction pipeline,
# collect all the relevant information, and will write it to a LSS catalogue
# file. 
#
# Current version has few columns and basically rewrites a zcat.fits file and a
# mtl.fits file but this will change as the code evolves.
#
# Lado Samushia March, 2017 (colkhis@gmail.com)

# The following information will be written into the LSS catalogue fits file:
# target_id -- unique target_id
# object_type -- this is one out of the following list: (GALAXY, QSO)
# redshift 
# ra -- right ascension
# dec -- declination
# redhsift_err -- redshift error
# redshift_warn -- redshift warning: 0 means ???, 4 means ???

from astropy.io import fits
from astropy.table import Table
import argparse
from datetime import date
import numpy as np
import os
from itertools import izip

# Parse arguments
# Two input files, one output file
# Input:
# zcat -- file containing measured redshifts
# mtl -- merged target list contains ra, dec
# targets -- file containing original targeting information
# (see http://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/mtl.html for
# the mtl file datamodel)
# Output:
# lsscat -- the output LSS catalogue
parser = argparse.ArgumentParser(description='Create a LSS catalogue.')
parser.add_argument('zcat',help='input redshift catalogue')
parser.add_argument('mtl',help='input merged target list file')
parser.add_argument('lsscat',help='output LSS catalogue')
args = parser.parse_args()

if os.path.isfile(args.lsscat) == True:
    print("This LSScat file already exists. Can not overwrite.")
    print("Choose an alternative name for the output LSScat file.")
    quit()

print("Extract information from mtl fits file")
mtl_file = fits.open(args.mtl)
mtl = mtl_file[1].data

mtl_target_id = mtl.field('targetid')
ra = mtl.field('ra')
dec = mtl.field('dec')
# desi_target 1 = LRG, 2 = ELG, 4 = QSO
desi_target = mtl.field('desi_target')
depth_r = mtl.field('depth_r')
galdepth_r = mtl.field('galdepth_r')
Ntarg = np.size(mtl_target_id)
# The arrays below will be filled in from the zcat file
z = np.zeros(Ntarg)
z_err = np.zeros(Ntarg)
z_warn = np.zeros(Ntarg)
stype = np.array(['000000' for item in range(Ntarg)])
# Create a dictionary from target_id to position in array
print("Create dictionary")
tid_dict = dict(zip(mtl_target_id,range(Ntarg)))

print("Extract information from zcat fits file")
zcat_file = fits.open(args.zcat)
zcat = zcat_file[1].data

zcat_target_id = zcat.field('targetid')
spec_type = zcat.field('spectype')
redshift = zcat.field('z')
redshift_err = zcat.field('zerr')
redshift_warn = zcat.field('zwarn')

print("select redshift and types")
for ID, st, zz, zz_err, zz_warn in izip(zcat_target_id, spec_type, redshift, redshift_err, redshift_warn):
    index = tid_dict[ID]
    stype[index] = st
    z[index] = zz
    z_err[index] = zz_err
    z_warn[index] = zz_warn

# Write information into LSS catalogue fits file
lss = Table([mtl_target_id, ra, dec, desi_target, depth_r, galdepth_r, stype, z, z_err, z_warn], names=('targetid','ra','dec','desitarget','depthr','galdepthr','stype','z','z_err','z_warn'))
lss.meta['comment'] = 'DESI Large-Scale Structure Catalogue'
lss.meta['date'] = date.today()
lss.write(args.lsscat, format='fits')
