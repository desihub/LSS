# This script will swipe through the output of the data reduction pipeline,
# collect all the relevant information, and will write it to a LSS catalogue
# file. 
#
# Current version has few columns and basically rewrites a zcat.fits file but
# this will change as the code evolves.
#
# Lado Samushia March, 2017 (colkhis@gmail.com)

# The following information will be written into the LSS catalogue fits file:
# target_id -- unique target_id
# object_type -- this is one out of the following list: (GALAXY, QSO)
# redshift 
# redhsift_err -- redshift error
# redshift_warn -- redshift warning: 0 means ???, 4 means ???

from astropy.io import fits
import argparse
import datetime

# Parse arguments
# Two input files, one output file
# Input:
# zcat -- file containing measured redshifts
# Output:
# lsscat -- the output LSS catalogue
parser = argparse.ArgumentParser(description='Create a LSS catalogue.')
parser.add_argument('zcat',help='input redshift catalogue')
parser.add_argument('lsscat',help='Output LSS catalogue')
args = parser.parse_args()

zcat_file = args.zcat
LSScat_file = args.lsscat

# Extract information from zcat fits file
hdulist = fits.open(zcat_file)
tbdata = hdulist[1].data
target_id = tbdata.field('targetid')

object_type = tbdata.field('spectype')
redshift = tbdata.field('z')
redshift_err = tbdata.field('zerr')
redshift_warn = tbdata.field('zwarn')

# Write information into LSS catalogue fits file
prihdr = fits.Header()
prihdr['comment'] = "DESI Large-Scale Structure Catalogue"
creation_date = datetime.today()
prihdr['creation date'] = creation_date.strftime('%m/%d/%Y')
prihdu = fits.PrimaryHDU(header=prihdr)

col1 = fits.Column(name='target_id',format='K',array=target_id)
col2 = fits.Column(name='object_type',format='10A',array=object_type)
col3 = fits.Column(name='redshift',format='D',array=redshift)
col4 = fits.Column(name='redshift_err',format='E',array=redshift_err)
col5 = fits.Column(name='redshift_warn',format='J',array=redshift_warn)

cols = fits.ColDefs([col1, col2, col3, col4, col5])
tnhdu = fits.BinTableHDU.from_coumns(cols)

thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto(LSScat_file)
