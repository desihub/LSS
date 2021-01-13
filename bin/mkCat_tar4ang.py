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

sys.path.append('../py')

#from this package
import LSS.imaging.select_samples as ss

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--tarver", help="version of targeting",default='0.44.0')
args = parser.parse_args()

type = args.type
tarver = args.tarver
version = '0' #integer for every tag that makes it to master

tp = 'DESI_TARGET'

outdir = '/project/projectdirs/desi/users/ajross/dr9/tarcat/v'+version+'/tv'+tarver+'/'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    print('created '+outdir)

dirsweeps = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/'
dirsweepn = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/'
targroot = '/project/projectdirs/desi/target/catalogs/dr9m/'+tarver+'/targets/main/resolve/'
ranroot =  '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-'
nran = 10

sfs = glob.glob(dirsweeps+'sweep*')
sfn = glob.glob(dirsweepn+'sweep*')



elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area; combination of bits applied to ELGs and LRGs in DR8 targeting

mkbsamp = False #make the base sample
domaskd = False #mask data based on mask bits above
domaskr = True #mask randoms

print('type being used for bright/dark '+type[:3])

if mkbsamp: #concatenate target files for given type, with column selection hardcoded
    prog = 'dark'
    if type[:3] == 'BGS':
        prog = 'bright'
    ss.gather_targets(type,targroot,outdir,tarver,prog)

if domaskd:
    dd = fitsio.read(outdir+type +'targetsDR9v'+tarver.strip('.')+'.fits'  )
    dd = ss.mask(dd,elgandlrgbits)
    outf = outdir+type +'targetsDR9v'+tarver.strip('.')+'_masked.fits'
    fitsio.write(outf,dd,clobber=True)
    print('wrote to '+outf)

if domaskr:     
    for ii in range(0,nran):
        rr = fitsio.read(ranroot+str(ii)+'.fits',columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS']) 
        #need to restrict columns on line above otherwise run out of memory
        rr = ss.mask(rr,elgandlrgbits)
        outf = outdir+'randomsDR9v'+tarver.strip('.')+'_'+str(ii)+'_masked.fits'
        fitsio.write(outf,rr,clobber=True)
        print('wrote to '+outf)
           

