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
targroot = '/project/projectdirs/desi/target/catalogs/dr9m/'+default+'/targets/main/resolve/'

sfs = glob.glob(dirsweeps+'sweep*')
sfn = glob.glob(dirsweepn+'sweep*')



elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area

mkbsamp = True #make the base sample

print('type being used for bright/dark '+type[:3])

if mkbsamp: #concatenate target files for given type, with column selection hardcoded
    prog = 'dark'
    if type[:3] == 'BGS':
        prog = 'bright'
    ss.gather_targets(type,targroot,outdir,tarver,prog)

