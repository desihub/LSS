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
#from desihub
#from desitarget import targetmask
#from regressis, must be installed
#from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--col_name", help="name of the column to add from data")
parser.add_argument("--replace", help="if the column is there, replace?",default='n')
parser.add_argument("--fix_weight", help="should replacing this column propagate to the WEIGHT column?",default='n')
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--blind", help="string to make output directory blinded or not",default='blinded/')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--inmode",help="string that denotes type of input file; '_clustering' would read from unsplit clustering catalog",default='_full_HPmapcut')
parser.add_argument("--par",help="process randoms in parallel?",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 


args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)

maindir = basedir +'/'+args.survey+'/LSS/'
ldirspec = maindir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'
dirin = dirout
dirout += args.blind

indata = Table(fitsio.read(dirin+args.tracer+args.inmode+'.dat.fits',columns=['TARGETID',args.col_name]))

regl = ['NGC','SGC']

for reg in regl:
    fname = dirout+args.tracer+'_'+reg+'_clustering.dat.fits'
    cd = Table(fitsio.read(fname))
    dojoin = 1
    if args.col_name in list(cd.dtype.names):
        if args.replace == 'y':
            if args.fix_weight == 'y':
                cd['WEIGHT'] /= cd[args.col_name]
            cd.remove_column(args.col_name)
        else:
            dojoin = 0
            #sys.exit('column is in catalog already! Set --replace y if you wish to replace it')
    if dojoin == 1:
        cd = join(cd,indata,keys=['TARGETID'],join_type='left')
        if args.fix_weight == 'y':
            cd['WEIGHT'] *= cd[args.col_name]
        common.write_LSS_scratchcp(cd,fname)
indata.rename_column('TARGETID', 'TARGETID_DATA')

def _add2ran(rn):
    for reg in regl:
        fname = dirout+args.tracer+'_'+reg+'_'+str(rn)+'_clustering.ran.fits'
        cd = Table(fitsio.read(fname))
        if args.col_name in list(cd.dtype.names):
            if args.replace == 'y':
                if args.fix_weight == 'y':
                    cd['WEIGHT'] /= cd[args.col_name]
                cd.remove_column(args.col_name)

            else:
                sys.exit('column is in random catalog '+str(rn)+' already! Use --replace y if you wish to replace it ')
        
        cd = join(cd,indata,keys=['TARGETID_DATA'],join_type='left')
        if args.fix_weight == 'y':
            cd['WEIGHT'] *= cd[args.col_name]

        common.write_LSS_scratchcp(cd,fname)

inds = np.arange(rm,rx)
if args.par == 'n':
    for rn in inds:
        _add2ran(rn)
if args.par == 'y':
	from multiprocessing import Pool
	with Pool() as pool:
		res = pool.map(_add2ran, inds)


#     for reg in regl:
#         fname = dirout+args.tracer+'_'+reg+'_'+str(rn)+'_clustering.ran.fits'
#         cd = Table(fitsio.read(fname))
#         if args.col_name in list(cd.dtype.names):
#             if args.replace == 'y':
#                 if args.fix_weight == 'y':
#                     cd['WEIGHT'] /= cd[args.col_name]
#                 cd.remove_column(args.col_name)
# 
#             else:
#                 sys.exit('column is in catalog already, but it was not in the data. Somthing strange happened! ')
#         
#         cd = join(cd,indata,keys=['TARGETID_DATA'],join_type='left')
#         if args.fix_weight == 'y':
#             cd['WEIGHT'] *= cd[args.col_name]
# 
#         common.write_LSS(cd,fname)
# 
