#script to read in existing LSS catalog information for some tracer and create a sub-sample
#Example submission to be added

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


import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main



parser = argparse.ArgumentParser()
parser.add_argument("--ccut",help="a string that is used define your subsample",default='FSFABSmagwecorr-R-20.5')
parser.add_argument("--input_tracer", help="tracer type that subsample will come from")
parser.add_argument("--basedir", help="base directory for input, default is SCRATCH",default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--outdir", help="directory for out, default is SCRATCH",default=os.environ['SCRATCH'])
parser.add_argument("--version", help="catalog version for input",default='v1.1')
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--use_map_veto", help="string to include in full file name denoting whether map veto was applied",default='_HPmapcut')
#parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='')

parser.add_argument("--compmd",help="use altmtl to use PROB_OBS for completeness weights in clustering catalogs",default='not_altmtl')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18,type=int) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--splitGC",help="convert to NGC/SGC catalogs",default='n')

parser.add_argument("--imsys",help="add weights for imaging systematics using eboss method?",default='n')
parser.add_argument("--imsys_clus",help="add weights for imaging systematics using eboss method, applied to clustering catalogs?",default='n')
parser.add_argument("--imsys_clus_ran",help="add weights for imaging systematics using eboss method, applied to clustering catalogs, to randoms?",default='n')
parser.add_argument("--imsys_clus_fb",help="perform linear weight fits in fine redshift bins",default='n')
parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=1,type=int)
parser.add_argument("--imsys_nside",help="healpix nside used for imaging systematic regressions",default=256,type=int)

parser.add_argument("--par", help="run different random number in parallel?",default='y')


args = parser.parse_args()
common.printlog(str(args),logger)

basedir = args.basedir
version = args.version
specrel = args.verspec
ccut = args.ccut
rm = int(args.minr)
rx = int(args.maxr)

#get some parameters that depend on the tracer
mainp = main(args.input_tracer,args.verspec,survey=args.survey) #parameters contained in mainp

dchi2 = mainp.dchi2 #used for good redshift selection when making clustering catalogs
#clustering catalogs will have these redshift bounds:
zmin = mainp.zmin
zmax = mainp.zmax



#set up input directory

maindir = basedir +'/'+args.survey+'/LSS/'

ldirspec = maindir+specrel+'/'

dirin = ldirspec+'LSScats/'+version+'/'
    
if not os.path.exists(dirin):
    
    sys.exit('issue with '+dirin+' it does not exist')

common.printlog('running subsampling '+dirin+args.input_tracer +' catalogs according to '+ccut,logger)

#parse arguments to toggle options
    
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

def get_FSF_loa(indata,fsf_cols,fsf_dir='/pscratch/sd/i/ioannis/fastspecfit/data/loa/catalogs/',prog='bright'):
    #add the fsf_cols to the existing data based on a TARGETID match
    #works with the data model that is new as of loa
    fsl = []
    for hp in range(0,12):
        fsi = fitsio.read(fsf_dir+'fastspec-loa-main-bright-nside1-hp'+str(hp).zfill(2)+'.fits',ext='SPECPHOT',columns = fsf_cols)
        fsl.append(fsi)
    fs = np.concatenate(fsl)
    del fsl
    ol = len(indata)
    indata = join(indata,fs,keys=['TARGETID']) #note, anything missing from fastspecfit will now be missing
    del fs
    common.printlog('length before/after fastspecfit join '+str(ol)+' '+str(len(indata)),logger)
    return indata

common.printlog('reading full data file '+dirin+args.input_tracer+'_full'+args.use_map_veto+'.dat.fits',logger)
fulldat = fitsio.read(dirin+args.input_tracer+'_full'+args.use_map_veto+'.dat.fits')

tracer_out = args.input_tracer+args.ccut

#for selections based on fastspecfit; other cases can be written similarly
if 'FSFABSmag' in args.ccut:
    #this is an example that can make subsamples based on fastspecfit absolute magnitudes
    #other critera can be added
    csplit = args.ccut.split('-')
    bnd = csplit[1]
    abmag = -float(csplit[2])
    fsf_cols = ['TARGETID','ABSMAG01_SDSS_'+bnd]
    #add more columns here based on args.ccut
    common.printlog('about to get columns from fastspecfit '+str(fsf_cols),logger)
    fulldat = get_FSF_loa(fulldat,fsf_cols)
    ecorr = np.zeros(len(fulldat))
    if 'ecorr' in args.ccut:
        ecorr = -0.8*(fulldat['Z_not4clus']-0.1) #seemed best here for getting constant n(z) /global/cfs/cdirs/desi/survey/catalogs/DA2/analysis/loa-v1/LSScats/BGS_explore.ipynb
    sel = fulldat['ABSMAG01_SDSS_'+bnd] < abmag + ecorr
    #add any additional selections here
    common.printlog('length after selection '+str(np.sum(sel)),logger)
    #write output to new "full" catalog at your defined location
    fout = args.outdir+'/'+tracer_out+'_full'+args.use_map_veto+'.dat.fits'
    common.write_LSS_scratchcp(fulldat[sel],fout,logger=logger)
    
    
#create "clustering" catalogs for data with no NGC/SGC split or FKP weights 
#needs to happen before randoms so randoms can get z and weights
weightileloc=True
if args.compmd == 'altmtl':
    weightileloc = False
if mkclusdat:
    ct.mkclusdat(args.outdir+'/'+tracer_out,weightileloc,tp=tracer_out,dchi2=dchi2,zmin=zmin,zmax=zmax,use_map_veto=args.use_map_veto)

#make clustering catalogs for randoms
nzcompmd = 'ran'
if args.compmd == 'altmtl':
    nzcompmd = args.compmd
rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL','TARGETID_DATA'] #columns to make sure are in the randoms
inds = np.arange(rm,rx)
if mkclusran:
    ranin = dirin + args.input_tracer 

    clus_arrays = [fitsio.read(args.outdir+tracer_out+'_clustering.dat.fits')]
    def _parfun_cr(ii):
        ct.mkclusran(ranin,dirout+'/'+tracer_out+'_',ii,rcols=rcols,clus_arrays=clus_arrays,use_map_veto=args.use_map_veto,compmd=nzcompmd,logger=logger,tp=args.input_tracer)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_parfun_cr, inds)

    else:
        for ii in inds:#range(rm,rx):
            _parfun_cr(ii)
        #,ntilecut=ntile,ccut=ccut)

