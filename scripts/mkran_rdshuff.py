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
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')



parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 


args = parser.parse_args()
print(args)

tp = args.tracer
basedir = args.basedir
version = args.version
specrel = args.verspec

rm = int(args.minr)
rx = int(args.maxr)


print('running catalogs for tracer type '+args.tracer)

    

if tp[:3] == 'BGS' or tp == 'bright' or tp == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.tracer,args.verspec,survey=args.survey)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tsnrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax



#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'


ldirspec = maindir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'
logfn = dirout+'log.txt'
if os.path.isfile(logfn):
    logf = open(logfn,'a')
else:
    logf = open(logfn,'w')
dirin = dirout

regl = ['_NGC','_SGC']

dz = .01
if tp[:3] == 'QSO':
	dz = 0.02
	zmin = 0.8
	zmax = 3.5
	P0 = 6000

if tp[:3] == 'LRG':
	P0 = 10000
	zmin = 0.4
	zmax = 1.1
if tp[:3] == 'ELG':
	P0 = 4000
	zmin = 0.8
	zmax = 1.6
if tp[:3] == 'BGS':
	P0 = 7000
	zmin = 0.1
	zmax = 0.4


for reg in regl:
    for ii in range(rm,rx):
        flin = dirout+tp+reg
        ffr = common.clusran_shufrd(flin,P0=P0,zmin=zmin,zmax=zmax,dz=dz)
        fnout = dirout+tp+'_rdshuf'+reg+'_'+str(ii)+'_clustering.ran.fits'
        common.write_LSS(ffr,fnout)