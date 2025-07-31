#Script to read in existing LSS catalog information for some tracer and create a sub-sample
#Expectation is that users will copy this and create their own criteria, following the example
#Example submission creates Mr < 20.5 catalog based on fastspecfit plus add hoc e correction with a further selection of 35% lowest star formation (percentile at which data looked bimodel)
#
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mkCat_subsamp.py --input_tracer BGS_ANY --mkfulldat y --clusd y --clusran y --nz y --splitGC y --ccut FSFABSmagwecorr-R-20.5-SFRlper-35 --imsys_clus y --imsys_clus_ran y
#Up to imaging systematics regression takes ~7 minutes ; imaging systematics takes another ~3 minutes
#
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

#LSS code, requires it git clone
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main



parser = argparse.ArgumentParser()
parser.add_argument("--ccut",help="a string that is used define your subsample",default='FSFABSmagwecorr-R-20.5-umzgper-50')
#arguments to find input data
parser.add_argument("--input_tracer", help="tracer type that subsample will come from")
parser.add_argument("--basedir", help="base directory for input, default is SCRATCH",default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--outdir", help="directory for out, default is SCRATCH",default=os.environ['SCRATCH'])
parser.add_argument("--version", help="catalog version for input",default='v1.1')
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--use_map_veto", help="string to include in full file name denoting whether map veto was applied",default='_HPmapcut')
#parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='')

parser.add_argument("--compmd",help="use altmtl to use PROB_OBS for completeness weights in clustering catalogs",default='not_altmtl')

#what steps to run (set all to y to get NGC/SGC clustering catalogs output)
parser.add_argument("--mkfulldat", help="whether to make the initial cut file that gets used throughout",default='n')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18,type=int) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--splitGC",help="convert to NGC/SGC catalogs",default='n')

#options for linear imaging systematic regressions
parser.add_argument("--imsys_clus",help="add weights for imaging systematics using eboss method, applied to clustering catalogs?",default='n')
parser.add_argument("--imsys_clus_ran",help="add weights for imaging systematics using eboss method, applied to clustering catalogs, to randoms?",default='n')
parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=1,type=int)
parser.add_argument("--usemaps", help="the list of maps to use; defaults to what is set by globals", type=str, nargs='*',default=None)
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

lssmapdirout = dirin+'/hpmaps/' #maps for imaging systematics regressions
    
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

def get_FSF_loa(indata,fsf_cols,fsf_dir='/dvs_ro/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/',prog='bright'):
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

tracer_out = args.input_tracer+args.ccut

if args.mkfulldat == 'y':

    common.printlog('reading full data file '+dirin+args.input_tracer+'_full'+args.use_map_veto+'.dat.fits',logger)
    fulldat = fitsio.read(dirin+args.input_tracer+'_full'+args.use_map_veto+'.dat.fits')
    
    
    
    #for selections based on fastspecfit; other cases can be written similarly
    if 'FSFABSmag' in args.ccut:
        #this is an example that can make subsamples based on fastspecfit absolute magnitudes
        #other critera can be added
        csplit = args.ccut.split('-')
        bnd = csplit[1]
        abmag = -float(csplit[2])
        fsf_cols = ['TARGETID','ABSMAG01_SDSS_'+bnd]
        #add more columns here based on args.ccut
        if 'umz' in args.ccut:
            fsf_cols.append('ABSMAG01_SDSS_U')
            fsf_cols.append('ABSMAG01_SDSS_Z')
            umz_str = csplit[3]
            umz_split = float(csplit[4]) #value to split on
            common.printlog('splitting on U-Z percentile '+str(umz_split),logger)
        if 'SFR' in args.ccut:
            fsf_cols.append('SFR')
            sfr_str = csplit[3]
            sfr_split = float(csplit[4]) #value to split on
            common.printlog('splitting on SFR percentile '+str(sfr_split),logger)
        common.printlog('about to get columns from fastspecfit '+str(fsf_cols),logger)
        fulldat = get_FSF_loa(fulldat,fsf_cols)
        ecorr = np.zeros(len(fulldat))
        if 'ecorr' in args.ccut:
            ecorr = -0.8*(fulldat['Z_not4clus']-0.1) #seemed best here for getting constant n(z) /global/cfs/cdirs/desi/survey/catalogs/DA2/analysis/loa-v1/LSScats/BGS_explore.ipynb
        sel = fulldat['ABSMAG01_SDSS_'+bnd] < abmag + ecorr
        common.printlog('length after Absmag selection '+str(np.sum(sel)),logger)
        if 'SFR' in args.ccut and 'per' in args.ccut: #'per' for percentile
            sel_sfr = fulldat['SFR'] > np.percentile(fulldat[sel]['SFR'],sfr_split)
            if 'g' in sfr_str: #'g' for greater than
                sel &= sel_sfr
            else:
                sel &= ~sel_sfr
            common.printlog('length after SFR selection '+str(np.sum(sel)),logger)
        if 'umz' in args.ccut and 'per' in args.ccut: #'per' for percentile
            sel_umz = (fulldat['ABSMAG01_SDSS_U']-fulldat['ABSMAG01_SDSS_Z']) > np.percentile((fulldat[sel]['ABSMAG01_SDSS_U']-fulldat[sel]['ABSMAG01_SDSS_Z']),umz_split)
            if 'g' in umz_str: #'g' for greater than
                sel &= sel_umz
            else:
                sel &= ~sel_umz
            common.printlog('length after UMZ selection '+str(np.sum(sel)),logger)

        #add any additional selections here
        
        #write output to new "full" catalog at your defined location
    else:
        sys.exit('should not have made it here, whatever you entered for --ccut did not trigger a cut, check code')
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
    ranin = dirin + args.input_tracer +'_'

    clus_arrays = [fitsio.read(args.outdir+'/'+tracer_out+'_clustering.dat.fits')]
    def _parfun_cr(ii):
        ct.mkclusran(ranin,args.outdir+'/'+tracer_out+'_',ii,rcols=rcols,clus_arrays=clus_arrays,use_map_veto=args.use_map_veto,compmd=nzcompmd,logger=logger,tp=args.input_tracer)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_parfun_cr, inds)

    else:
        for ii in inds:#range(rm,rx):
            _parfun_cr(ii)

#define P0 value used for fiducial FKP weights and dz used for creating n(z)
if tracer_out[:3] == 'QSO':
    dz = 0.02
    P0 = 6000    
else:    
    dz = 0.01

if tracer_out[:3] == 'LRG':
    P0 = 10000
if tracer_out[:3] == 'ELG':
    P0 = 4000
if tracer_out[:3] == 'BGS':
    P0 = 7000

nran = rx-rm
regions = ['NGC', 'SGC']

#function to take a file and split it NGC/SGC
def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    fn = Table(fitsio.read(flroot.replace('global','dvs_ro') +app))
    sel_ngc = common.splitGC(fn)#gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS_scratchcp(fn[sel_ngc],outf_ngc,logger=logger)
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS_scratchcp(fn[~sel_ngc],outf_sgc,logger=logger)


dirout = args.outdir #just because original copied code used dirout

#split catalogs NGC/SGC
if args.splitGC == 'y':
    fb = dirout+'/'+tracer_out+'_'
    splitGC(fb,'.dat')
    def _spran(rann):
        splitGC(fb,'.ran',rann)
    inds = np.arange(nran)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_spran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _spran(rn)

#this 1) calculates the n(z) for the sample
#2) adds the n(z) as a function of completeness to the catalogs (column NX)
#3) calculates the FKP weights based on NX
#4) refactors the weights (section 8.2 of the KP3 paper arXiv:2411.12020)
#5) writes the NGC/SGC clustering catalogs back out for data and randoms
#It needs to take inputs for completeness and mean weight as a function of NTILE, of the non-subsampled catalog, so that the angular upweighting option can remain consistent
def get_ntile_info(fd):
    ntl = np.unique(fd['NTILE'])
    comp_ntl = np.ones(len(ntl))
    weight_ntl = np.ones(len(ntl))
    for i in range(0,len(ntl)):
        sel = fd['NTILE'] == ntl[i]
        mean_ntweight = np.mean(fd['WEIGHT_COMP'][sel])        
        weight_ntl[i] = mean_ntweight
        comp_ntl[i] = 1/mean_ntweight#*mean_fracobs_tiles
        
        if compmd != 'altmtl':
            fttl = np.zeros(len(ntl))
            for i in range(0,len(ntl)): 
                sel = fd['NTILE'] == ntl[i]
                mean_fracobs_tiles = np.mean(fd[sel]['FRAC_TLOBS_TILES'])
                fttl[i] = mean_fracobs_tiles
        else:
            fttl = np.ones(len(ntl))
    comp_ntl = comp_ntl*fttl
    return comp_ntl,weight_ntl
 
if args.nz == 'y':
    for reg in regions:#allreg:
        #file names
        fb = dirout+'/'+tracer_out+'_'+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.txt'
        #make n(z)
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,compmd=nzcompmd)
        #do steps 2-5 above
        extra_dir = 'nonKP'
        if args.compmd == 'altmtl':
            extra_dir = 'PIP'
        clus_orig = fitsio.read(dirin+'/'+extra_dir+'/'+args.input_tracer+'_'+reg+'_clustering.dat.fits')
        comp_ntl,weight_ntl = get_ntile_info(clus_orig)                        
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,nran=nran,par=args.par,compmd=nzcompmd,comp_ntl=comp_ntl,weight_ntl=weight_ntl,logger=logger)

# determine linear weights for imaging systematics
# this is new for doing after the fact based on clustering catalogs
if args.imsys_clus == 'y':
    #import package
    from LSS.imaging import densvar 
    
    #setup redshift bins to use for regressions; you might find something else works better!
    if args.input_tracer[:3] == 'ELG':
        if args.imsys_zbin == 'y':
            zrl = [(0.8,1.1),(1.1,1.6)]
        else:
            zrl = [(0.8,1.6)]
    if args.input_tracer[:3] == 'QSO':
        if args.imsys_zbin == 'y':
            zrl = [(0.8,1.3),(1.3,2.1),(2.1,3.5)] 
        else:
            zrl = [(0.8,3.5)]   
    if args.input_tracer[:3] == 'LRG':
        if args.imsys_zbin == 'y':
            zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
        else:
            zrl = [(0.4,1.1)]  
    if args.input_tracer[:3] == 'BGS':
        zrl = [(0.1,0.4)]

    
    #get maps to regress against
    if args.usemaps == None:
        fit_maps = mainp.fit_maps
    else:
        fit_maps = [mapn for mapn in args.usemaps]

    #get NGC/SGC catalogs and stack them (weights will be fit splitting the data into different photometric regions)
    fname = os.path.join(dirout, tracer_out+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname))
    fname = os.path.join(dirout, tracer_out+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname))
    dat = vstack([dat_sgc,dat_ngc])
    #foutname = os.path.join(dirout, tracer_clus+'_clustering.dat.fits')
    #get randoms
    ranl = []
    for i in range(0,args.nran4imsys):
        ran = fitsio.read(os.path.join(dirout, tracer_out+'_NGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)
        ran = fitsio.read(os.path.join(dirout, tracer_out+'_SGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)
    rands = np.concatenate(ranl)
    
    #fiducial column name for weights, initialize as 1.
    syscol = 'WEIGHT_IMLIN_CLUS' 
    dat[syscol] = np.ones(len(dat))
    #photometric regions
    regl = ['S','N']
    if args.input_tracer == 'QSO':
        regl = ['DES','SnotDES','N']
    
    #do fit looping over regions and redshift bins
    for reg in regl:
        #handling for loading maps of potential systematics to regress against
        regu = reg
        if reg == 'DES' or reg == 'SnotDES':
            regu = 'S'
        mptr = args.input_tracer
        if args.input_tracer == 'BGS_ANY':
            mptr = 'BGS_BRIGHT' #bright and any cover the same footprint so the same map is used for both
        pwf = lssmapdirout+mptr+'_mapprops_healpix_nested_nside256_'+regu+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            #apply extinction corrections to depth to get total depth estimate
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        #Delta EBV using DESI EBV maps
        debv = common.get_debv()
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        
        selr = rands['PHOTSYS'] == reg

        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            
            common.printlog('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax),logger)
            if args.input_tracer == 'LRG' and args.usemaps == None:
                if reg == 'N':
                    fitmapsbin = fit_maps
                else:
                    if zmax == 0.6:
                        fitmapsbin = mainp.fit_maps46s
                    if zmax == 0.8:
                        fitmapsbin = mainp.fit_maps68s
                    if zmax == 1.1:
                        fitmapsbin = mainp.fit_maps81s
            else:
                fitmapsbin = fit_maps
            use_maps = fitmapsbin
            #now, everything in place to actually perform regression for this reg and zbin
            wsysl = densvar.get_imweight(dat,rands,zmin,zmax,reg,fitmapsbin,use_maps,sys_tab=sys_tab,zcol='Z',modoutname = dirout+'/'+tracer_out+'_'+reg+'_'+str(zmin)+str(zmax)+'_linfitparam.txt',figname=dirout+'/'+tracer_out+'_'+reg+'_'+str(zmin)+str(zmax)+'_linclusimsysfit.png',wtmd='clus')
            #we want to update the weights for the selection of data just input to the regression
            sel = wsysl != 1 
            dat[syscol][sel] = wsysl[sel]
    #attach data to NGC/SGC catalogs, write those out
    #we will do a join
    dat.keep_columns(['TARGETID',syscol])
    
    if syscol in dat_ngc.colnames:
        dat_ngc.remove_column(syscol)
    dat_ngc = join(dat_ngc,dat,keys=['TARGETID'])
    #apply weight to final weight columns, remove any previous weighting
    dat_ngc['WEIGHT'] /= dat_ngc['WEIGHT_SYS']
    dat_ngc['WEIGHT_SYS'] = dat_ngc[syscol]
    dat_ngc['WEIGHT'] *= dat_ngc['WEIGHT_SYS']
    #write out NGC
    common.write_LSS_scratchcp(dat_ngc,os.path.join(dirout, tracer_out+'_NGC_clustering.dat.fits'),logger=logger)
    #do SGC
    if syscol in dat_sgc.colnames:
        dat_sgc.remove_column(syscol)
    dat_sgc = join(dat_sgc,dat,keys=['TARGETID'])
    #apply weight to final weight columns
    dat_sgc['WEIGHT'] /= dat_sgc['WEIGHT_SYS']
    dat_sgc['WEIGHT_SYS'] = dat_sgc[syscol]
    dat_sgc['WEIGHT'] *= dat_sgc['WEIGHT_SYS']
    #write out SGC
    common.write_LSS_scratchcp(dat_sgc,os.path.join(dirout, tracer_out+'_SGC_clustering.dat.fits'),logger=logger)

#column needs to be added to randoms
if args.imsys_clus_ran == 'y':
    #do randoms
    syscol = 'WEIGHT_IMLIN_CLUS'
    fname = os.path.join(dirout, tracer_out+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    fname = os.path.join(dirout, tracer_out+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    dat = vstack([dat_sgc,dat_ngc])
    dat.rename_column('TARGETID','TARGETID_DATA') #randoms have their weights modulated based on the data used for the redshift
    regl = ['NGC','SGC']
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(dirout, tracer_out+'_'+reg+'_'+str(rann)+'_clustering.ran.fits')
            ran = Table(fitsio.read(ran_fn))
            if syscol in ran.colnames:
                ran.remove_column(syscol)
            ran = join(ran,dat,keys=['TARGETID_DATA'])
            ran['WEIGHT'] /= ran['WEIGHT_SYS'] #remove effect of any original weighting
            ran['WEIGHT'] *= ran[syscol]
            ran['WEIGHT_SYS'] = ran[syscol]
            common.write_LSS_scratchcp(ran,ran_fn,logger=logger)

    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _add2ran(rn)

