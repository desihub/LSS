'''
one executable to create catalogs for given target type meant for angular clustering
'''



#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt

#sys.path.append('../py')

#from this package
import LSS.imaging.select_samples as ss

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--tarver", help="version of targeting",default='0.57.0')
parser.add_argument("--survey", help="e.g., sv1 or main",default='sv3')
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()

type = args.type
tarver = args.tarver
version = args.version
basedir = args.basedir
survey = args.survey

if survey == 'main':
    tp = 'DESI_TARGET'
    sw = ''
if survey == 'sv1':
    tp = 'SV1_DESI_TARGET'
    sw = 'sv1'
if survey == 'sv3':
    tp = 'SV3_DESI_TARGET'
    sw = 'sv3'

outdir = basedir+'/tarcat/v'+version+'/tv'+tarver+'/'
if not os.path.exists( basedir+'/tarcat'):
    os.mkdir(basedir+'/tarcat')
    print('created '+basedir+'/tarcat')

if not os.path.exists( basedir+'/tarcat/v'+version):
    os.mkdir(basedir+'/tarcat/v'+version)
    print('created '+basedir+'/tarcat/v'+version)

if not os.path.exists(outdir):
    os.mkdir(outdir)
    print('created '+outdir)

dirsweeps = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/'
dirsweepn = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/'
targroot = '/project/projectdirs/desi/target/catalogs/dr9/'+tarver+'/targets/'+survey+'/resolve/'
ranroot =  '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-'
nran = 10

sfs = glob.glob(dirsweeps+'sweep*')
sfn = glob.glob(dirsweepn+'sweep*')



elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area; combination of bits applied to ELGs and LRGs in DR8 targeting

mkbsamp = True #make the base sample
domaskd = True #mask data based on mask bits above
domaskr = True #mask randoms
'test'
print('type being used for bright/dark '+type[:3])

#columns to select from target sample
keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
	   'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
	   'MASKBITS', 'EBV', 'PHOTSYS','TARGETID',tp,'SHAPE_R']


if mkbsamp: #concatenate target files for given type, with column selection hardcoded
    prog = 'dark'
    if type[:3] == 'BGS':
        prog = 'bright'
    ss.gather_targets(type,targroot,outdir,tarver,survey,prog,keys=keys)

if domaskd:
    dd = fitsio.read(outdir+type+sw +'targetsDR9v'+tarver.strip('.')+'.fits'  )
    dd = ss.mask(dd,elgandlrgbits)
    outf = outdir+type+sw +'targetsDR9v'+tarver.strip('.')+'_masked.fits'
    fitsio.write(outf,dd,clobber=True)
    print('wrote to '+outf)

if domaskr:     
    for ii in range(0,nran):
        rr = fitsio.read(ranroot+str(ii)+'.fits',columns=['RA','DEC','BRICKID','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS']) 
        #need to restrict columns on line above otherwise run out of memory
        rr = ss.mask(rr,elgandlrgbits)
        outf = outdir+'randomsDR9v'+tarver.strip('.')+'_'+str(ii)+'_masked.fits'
        fitsio.write(outf,rr,clobber=True)
        print('wrote to '+outf)
           

