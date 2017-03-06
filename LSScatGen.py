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
import argparse
from datetime import date
import numpy as np

# Parse arguments
# Two input files, one output file
# Input:
# zcat -- file containing measured redshifts
# mtl -- merged target list contains ra, dec
# (see http://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/mtl.html for
# the mtl file datamodel)
# Output:
# lsscat -- the output LSS catalogue
parser = argparse.ArgumentParser(description='Create a LSS catalogue.')
parser.add_argument('zcat',help='input redshift catalogue')
parser.add_argument('mtl',help='input merged target list file')
parser.add_argument('lsscat',help='output LSS catalogue')
args = parser.parse_args()

zcat_file = args.zcat
mtl_file = args.mtl
LSScat_file = args.lsscat

print("Extract information from zcat fits file")
hdulist = fits.open(zcat_file)
tbdata = hdulist[1].data

target_id = tbdata.field('targetid')
object_type = tbdata.field('spectype')
redshift = tbdata.field('z')
redshift_err = tbdata.field('zerr')
redshift_warn = tbdata.field('zwarn')

hdulist.close()

print("Apply selection")
# Binary mask to only select GALAXY or QSO
isqso = object_type == 'QSO'
isgalaxy = object_type == 'GALAXY'
bm = np.any([isqso,isgalaxy],axis=0)
target_id = target_id[bm]
object_type = object_type[bm]
redshift = redshift[bm]
redshift_err = redshift_err[bm]
redshift_warn = redshift_warn[bm]

print("Extract information from mtl fits file")
hdulist = fits.open(mtl_file)
tbdata = hdulist[1].data

# The ordering of targets in the mtl file is not necessarily the same as the
# ordering in the zcat file. Also, we have applied some masks.
mtl_target_id = tbdata.field('targetid')
mtl_ra = tbdata.field('ra')
mtl_dec = tbdata.field('dec')

# Create python dictionary to select ra, dec for a specific target_id
print("Create radec dictionary")
radec = zip(mtl_ra, mtl_dec)
radec_dict = dict(zip(mtl_target_id, radec))

ra = np.zeros(np.size(target_id))
dec = np.zeros(np.size(target_id))

print("select ra, dec")
for index, tid in enumerate(target_id):
    ra[index], dec[index] = radec_dict[tid]

# Write information into LSS catalogue fits file
print("write to LSS catalogue file")
prihdr = fits.Header()
prihdr['comment'] = "DESI Large-Scale Structure Catalogue"
creation_date = date.today()
prihdr['date'] = creation_date.strftime('%m/%d/%Y')
prihdu = fits.PrimaryHDU(header=prihdr)

col1 = fits.Column(name='target_id',format='K',array=target_id)
col2 = fits.Column(name='object_type',format='10A',array=object_type)
col3 = fits.Column(name='redshift',format='D',array=redshift)
col4 = fits.Column(name='ra',format='D',array=ra)
col5 = fits.Column(name='dec',format='D',array=dec)
col6 = fits.Column(name='redshift_err',format='E',array=redshift_err)
col7 = fits.Column(name='redshift_warn',format='J',array=redshift_warn)

cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7])
tbhdu = fits.BinTableHDU.from_columns(cols)

thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto(LSScat_file)
