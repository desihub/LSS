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

import logging
logname = 'mkCat'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


#logger = logging.getLogger('mkCat')
#logger.setLevel(level=logging.INFO)

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
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--basedir_blind", help="base directory for output for blinded catalogs, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--redotar", help="remake the target file for the particular type (needed if, e.g., the requested columns are changed)",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='n')
parser.add_argument("--add_veto", help="add veto column for given type, matching to targets",default='n')
parser.add_argument("--join_etar", help="whether or not to join to the target files with extra brick pixel info",default='n')
parser.add_argument("--apply_veto", help="apply vetos for imaging, priorities, and hardware failures",default='n')
parser.add_argument("--mkHPmaps", help="make healpix maps for imaging properties using sample randoms",default='n')
parser.add_argument("--usemaps", help="the list of maps to use; defaults to what is set by globals", type=str, nargs='*',default=None)
parser.add_argument("--apply_map_veto", help="apply vetos to data and randoms based on values in healpix maps",default='n')
parser.add_argument("--use_map_veto", help="string to include in full file name denoting whether map veto was applied",default='_HPmapcut')
parser.add_argument("--add_tlcomp", help="add completeness FRAC_TLOBS_TILES to randoms",default='n')



parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='')

parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18,type=int) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--nzfull", help="get n(z) from full files",default='n')

parser.add_argument("--FKPfull", help="add FKP weights to full catalogs",default='n')
parser.add_argument("--addnbar_ran", help="just add nbar/fkp to randoms",default='n')
parser.add_argument("--add_ke", help="add k+e corrections for BGS data to clustering catalogs",default='n')
parser.add_argument("--add_fs", help="add rest frame info from fastspecfit",default='n')
parser.add_argument("--redoBGS215", help="whether to use purely photometry+z based rest frame info or fastspecfit",default='n')
parser.add_argument("--absmagmd", help="whether to use purely photometry+z based rest frame info or fastspecfit, or no k correction at all",choices=['spec','phot','nok'],default='spec')

parser.add_argument("--blinded", help="are we running on the blinded full catalogs?",default='n')

parser.add_argument("--swap20211212", help="swap petal 9 redshifts from 20211212",default='n')
parser.add_argument("--fixzwarn", help="change any originally 'not observed' zwarn back to 999999",default='n')


parser.add_argument("--prepsysnet",help="prepare data to get sysnet weights for imaging systematics?",default='n')
parser.add_argument("--add_sysnet",help="add sysnet weights for imaging systematics to full files?",default='n')
parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='y')

parser.add_argument("--imsys",help="add weights for imaging systematics using eboss method?",default='n')
parser.add_argument("--imsys_clus",help="add weights for imaging systematics using eboss method, applied to clustering catalogs?",default='n')
parser.add_argument("--imsys_clus_ran",help="add weights for imaging systematics using eboss method, applied to clustering catalogs, to randoms?",default='n')
parser.add_argument("--imsys_clus_fb",help="perform linear weight fits in fine redshift bins",default='n')



parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=1,type=int)

parser.add_argument("--regressis",help="RF weights for imaging systematics?",default='n')
parser.add_argument("--regmode",help="RF and Linear are choices",default='RF')
parser.add_argument("--add_regressis",help="add RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis_ext",help="add RF weights for imaging systematics, calculated elsewhere",default='n')
parser.add_argument("--imsys_nside",help="healpix nside used for imaging systematic regressions",default=256,type=int)
parser.add_argument("--imsys_colname",help="column name for fiducial imaging systematics weight, if there is one (array of ones by default)",default=None)


parser.add_argument("--add_weight_zfail",help="add weights for redshift systematics to full file?",default='n')
parser.add_argument("--add_bitweight",help="add info from the alt mtl",default='n')
parser.add_argument("--compmd",help="use altmtl to use PROB_OBS",default='not_altmtl')
parser.add_argument("--addNtileweight2full",help="whether to add the NTILE weight to the full catalogs (necessary for consistent angular upweighting)",default='n')
parser.add_argument("--NStoGC",help="convert to NGC/SGC catalogs",default='n')
parser.add_argument("--splitGC",help="convert to NGC/SGC catalogs",default='n')
parser.add_argument("--resamp",help="resample radial info for different selection function regions",default='n')


parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

parser.add_argument("--par", help="run different random number in parallel?",default='y')

parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)
parser.add_argument("--ranonly",help="if y, only operate on randoms when applying vetos",default='n')
parser.add_argument("--readpars",help="set to y for certain steps if you want to read from previous fits",default='n')


#options not typically wanted

parser.add_argument("--ran_utlid", help="cut randoms so that they only have 1 entry per tilelocid",default='n')
parser.add_argument("--swapz", help="if blinded, swap some fraction of redshifts?",default='n')

#AJRM 
parser.add_argument("--use_allsky_rands", help="if yes, use all sky randoms to get fractional area per pixel for SYSNet data preparation",default='n')

args = parser.parse_args()
common.printlog(str(args),logger)

type = args.type
tp = type
basedir = args.basedir
version = args.version
specrel = args.verspec
ntile = args.ntile
ccut = args.ccut
rm = int(args.minr)
rx = int(args.maxr)


common.printlog('running catalogs for tracer type '+type,logger)

redotar = False
if args.redotar == 'y':
    redotar = True

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
        
    
if mkfulld:
    common.printlog('making "full" catalog file for data',logger)    
    
    
#mkfullr = True #make the random files associated with the full data files
#if args.fullr == 'n':
#    mkfullr = False
    
#if mkfullr:
#    print('making full catalog for randoms, files '+str(rm)+ ' through '+str(rx))
#    print('(if running all, consider doing in parallel)')    
    
mkclusdat = False
mkclusran = False
if args.clusd == 'y':
    mkclusdat = True
    
if mkclusdat:
    common.printlog('making clustering catalog for data',logger)
    
if args.clusran == 'y':
    mkclusran = True
    
if mkclusran:
    common.printlog('making clustering catalog for randoms, files '+str(rm)+ ' through '+str(rx),logger)
    common.printlog('(if running all, consider doing in parallel)',logger)  

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.type,args.verspec,survey=args.survey)
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


#wt = mtld['FAPRGRM'] == progl
#wt &= mtld['SURVEY'] == 'main'
#wt &= mtld['ZDONE'] == 'true'
#mtld = mtld[wt]
#print('there are '+str(len(mtld))+' tiles')


#columns to select from target sample
keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
       'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
       'MASKBITS','WISEMASK_W1','WISEMASK_W2', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET','BGS_TARGET','SHAPE_R']

if type[:3] == 'MWS':
    keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
       'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
       'MASKBITS','WISEMASK_W1','WISEMASK_W2', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET','MWS_TARGET','SHAPE_R']


#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

#if not os.path.exists(maindir+'/logs'):
#    os.mkdir(maindir+'/logs')
#    print('made '+maindir+'/logs')

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    common.printlog('made '+ldirspec,logger)
    
if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    common.printlog('made '+ldirspec+'LSScats',logger)

dirout = ldirspec+'LSScats/'+version+'/'

if not os.path.exists(dirout):
    os.mkdir(dirout)
    common.printlog('made '+dirout,logger)

if not os.path.exists(dirout+args.extra_clus_dir):
    os.mkdir(dirout+args.extra_clus_dir)
    common.printlog('made '+dirout+args.extra_clus_dir,logger)


args.extra_clus_dir
logfn = dirout+'log.txt'
if os.path.isfile(logfn):
    logf = open(logfn,'a')
else:
    logf = open(logfn,'w')
dirin = dirout
if args.blinded == 'y':
    dirout = args.basedir_blind+'/LSScats/'+version+'/blinded/'
    #dirout += 'blinded/'

if not os.path.exists(dirout):
    os.mkdir(dirout)
    logger.printlog('made '+dirout,logger)    

tarver = '1.1.1'
tardir = '/global/cfs/cdirs/desi/target/catalogs/dr9/'+tarver+'/targets/main/resolve/'
tarf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+type +'targetsDR9v'+tarver.strip('.')+'.fits'

mktar = True
if os.path.isfile(tarf) and redotar == False or len(type.split('-'))>1:
    mktar = False
#if type == 'BGS_BRIGHT':
#    mktar = False    
