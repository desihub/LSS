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
parser.add_argument("--tracer", help="tracer type to be selected",default='ELG_LOPnotqso')
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--maskbits", help="string identifying the maskbits to use",default='all')
parser.add_argument("--blinded", help="string for blinded directory",default='blinded/')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--par", help="use parallel processing", default='y')

args = parser.parse_args()
print(args)

tracer = args.tracer
tp = tracer
basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)


print('running catalogs for tracer type '+tp)

    

    

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

#if not os.path.exists(maindir+'/logs'):
#    os.mkdir(maindir+'/logs')
#    print('made '+maindir+'/logs')

ldirspec = maindir+specrel+'/'
    

dirout = ldirspec+'LSScats/'+version+'/'
logfn = dirout+'log.txt'
if os.path.isfile(logfn):
    logf = open(logfn,'a')
else:
    logf = open(logfn,'w')
dirin = dirout



tracer_clus = tracer #legacy of previous script

datin = fitsio.read(dirout + tracer+'_full_noveto.dat.fits',columns=['TARGETID','elg_mask'])
if args.maskbits == 'all':
    inmask = datin['elg_mask'] > 0
    mask_tids = np.unique(datin[inmask]['TARGETID'])
print('got masked targetids')
regl = ['NGC','SGC']
for reg in regl:
    fn = dirout + args.blinded + tracer +'_'+reg+'_clustering.dat.fits'
    fno = dirout + args.blinded + tracer+'elgmask_'+args.maskbits +'_'+reg+'_clustering.dat.fits'
    clusin = fitsio.read(fn)
    masked = np.isin(clusin['TARGETID'],mask_tids)
    print(reg,len(clusin),len(clusin[~masked]))
    common.write_LSS(clusin[~masked],fno)

def _parfun(rannum):
    datin = fitsio.read(dirout + tracer+'_'+str(rannum)+'_full_noveto.ran.fits',columns=['TARGETID','elg_mask'])
    if args.maskbits == 'all':
        inmask = datin['elg_mask'] > 0
        mask_tids = np.unique(datin[inmask]['TARGETID'])
    print(rannum,'got masked targetids')
    regl = ['NGC','SGC']
    for reg in regl:
        fn = dirout + args.blinded + tracer +'_'+reg+'_'+str(rannum)+'_clustering.ran.fits'
        fno = dirout + args.blinded + tracer+'elgmask_'+args.maskbits +'_'+reg+'_'+str(rannum)+'_clustering.ran.fits'
        clusin = fitsio.read(fn)
        masked = np.isin(clusin['TARGETID'],mask_tids)
        print(reg,len(clusin),len(clusin[~masked]))
        common.write_LSS(clusin[~masked],fno)

nran = args.maxr-args.minr
inds = np.arange(args.minr,args.maxr)
if args.par == 'y':    
    from multiprocessing import Pool
    with Pool(processes=nran*2) as pool:
        res = pool.map(_parfun, inds)
else:
    for ii in inds:
        _parfun(ii)    


