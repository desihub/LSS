'''
This script is meant to take a mock with coordinates RA,DEC,Z (but with arbitrary names in inputs)
and convert it to occupy the footprint of some set of tiles (one of the input variables)
and match the data model of LSS catalogs

example run
srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python scripts/mock_tools/skymock2LSSclus_fast.py --realization 0 
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
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
from LSS.globals import main

import logging
logger = logging.getLogger('mkCat')
logger.setLevel(level=logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

logger.info('run started')

#import LSS.mkCat_singletile.fa4lsscat as fa
#from LSS.globals import main

if os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    logger.info('NERSC_HOST is not permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
#parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--realization",type=int)
parser.add_argument("--prog", default="DARK")
parser.add_argument("--survey",default='DA2')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is all 18)",default=18,type=int) 
parser.add_argument("--mockver", default='AbY3HF', help = "which mocks to use")

parser.add_argument("--tracer", default = 'LRG')
parser.add_argument("--snapshot", default = 'z0.725')
parser.add_argument("--outloc",help='default will write to your scratch', default = None)
parser.add_argument("--par", default = 'y',help='whether to run random steps in parallel or not')
parser.add_argument("--mkdat", default = 'y')
parser.add_argument("--mkran", default = 'y')
parser.add_argument("--nz", default = 'y')
parser.add_argument("--splitGC", default = 'y')
#parser.add_argument("--remove_unassigned", default = 'y', help = 'set to n if dont want to include unassigned targets in catalog')


args = parser.parse_args()
logger.info(args)
tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/tiles-{PROG}.fits'.format(PROG = args.prog.upper()))
rm = int(args.minr)
rx = int(args.maxr)

notqso = ''

if args.mockver == 'AbY3HF':
    base_dir = '/dvs_ro/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph0'
    mockdir = base_dir+str(args.realization).zfill(2)+'/CutSky/'+args.tracer+'/'+args.snapshot+'/'
    if args.mkdat == 'y':
        in_data_fn = mockdir+'/'+'cutsky_'+args.tracer+'_'+args.snapshot+'_AbacusSummit_base_c000_ph0'+str(args.realization).zfill(2)+'.fits'
        in_data_fn = in_data_fn.replace('global','dvs_ro')
        logger.info(in_data_fn)
        cols = ['RA','DEC','Z']
        mock_data = Table(fitsio.read(in_data_fn,columns=cols))
        nin = len(mock_data)
        mock_data['TARGETID'] = np.arange(nin).astype(int)
    
        nuid = len(np.unique(mock_data['TARGETID'])) #check that we really have unique targetid; I think a type of int should always be large enough
        if nuid != nin:
            sys.exit('TARGETID are not unique!')
        selfoot = is_point_in_desi(tiletab,mock_data['RA'],mock_data['DEC'])
        mock_data = mock_data[selfoot]
        logger.info('length before/after cut to footprint '+str(nin)+'/'+str(len(mock_data)))
        mock_data = common.addNS(mock_data)
        logger.info('numbers in different photometric regions '+str(np.unique(mock_data['PHOTSYS'],return_counts=True)))
    tracerd = args.tracer+args.snapshot

if args.mockver == 'Uchuu-SHAM_Y3':
    base_dir = '/dvs_ro/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/complete/'
    mockdir = base_dir
    if args.mkdat == 'y':
        in_data_fn = mockdir+'/'+'Uchuu-SHAM_'+args.tracer+'_Y3-v2.0_'+str(args.realization).zfill(4)+'_clustering.dat.fits'
        in_data_fn = in_data_fn.replace('global','dvs_ro')
        logger.info(in_data_fn)
        cols = ['RA','DEC','Z']
        mock_data = Table(fitsio.read(in_data_fn,columns=cols))
        nin = len(mock_data)
        mock_data['TARGETID'] = np.arange(nin).astype(int)
    
        nuid = len(np.unique(mock_data['TARGETID'])) #check that we really have unique targetid; I think a type of int should always be large enough
        if nuid != nin:
            sys.exit('TARGETID are not unique!')
        selfoot = is_point_in_desi(tiletab,mock_data['RA'],mock_data['DEC'])
        mock_data = mock_data[selfoot]
        logger.info('length before/after cut to footprint '+str(nin)+'/'+str(len(mock_data)))
        mock_data = common.addNS(mock_data)
        logger.info('numbers in different photometric regions '+str(np.unique(mock_data['PHOTSYS'],return_counts=True)))
    tracerd = args.tracer

if args.mockver == 'holiv2':
    base_dir = '/global/cfs/cdirs/desi/mocks/cai/holi/v2.0/'
    mockdir = base_dir +'/seed'+str(args.realization).zfill(4)
    if args.mkdat == 'y':
        import h5py
        import hdf5plugin #need to be in the cosmodesi test environment, as of Sep 4th 25
        in_data_fn = mockdir+'/holi_'+args.tracer+'_v2.0_GCcomb_clustering.dat.h5'
        mock_data = Table()
        with h5py.File(in_data_fn) as fn:
            
            columns = fn.keys()
            for col in columns:
                mock_data[col.upper()] = fn[col][:]
        nin = len(mock_data)
        mock_data['TARGETID'] = np.arange(nin).astype(int)
    
        nuid = len(np.unique(mock_data['TARGETID'])) #check that we really have unique targetid; I think a type of int should always be large enough
        if nuid != nin:
            sys.exit('TARGETID are not unique!')
        selfoot = is_point_in_desi(tiletab,mock_data['RA'],mock_data['DEC'])
        mock_data = mock_data[selfoot]
        logger.info('length before/after cut to footprint '+str(nin)+'/'+str(len(mock_data)))
        mock_data = common.addNS(mock_data)
        logger.info('numbers in different photometric regions '+str(np.unique(mock_data['PHOTSYS'],return_counts=True)))
    tracerd = args.tracer



if args.outloc == None:
    outdir = os.getenv(scratch)+'/'+args.mockver+'/mock'+str(args.realization)+'/'

else: 
    outdir = args.outloc+'/'+args.mockver+'/mock'+str(args.realization)+'/'


if not os.path.exists(outdir):
    os.makedirs(outdir)


logger.info('input directory is '+mockdir)
logger.info('output directory is '+outdir)    


out_data_fn = outdir+tracerd+'_tilearea_clustering.dat.fits'
out_data_froot = outdir+tracerd+'_tilearea_'



def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    fn = Table(fitsio.read(flroot.replace('global','dvs_ro') +app))
    #if datran == '.ran':
    #    fn.keep_columns(['RA', 'DEC', 'Z', 'WEIGHT', 'WEIGHT_FKP', 'TARGETID_DATA'])
    #c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    #gc = c.transform_to('galactic')
    sel_ngc = common.splitGC(fn)#gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS_scratchcp(fn[sel_ngc],outf_ngc,logger=logger)
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS_scratchcp(fn[~sel_ngc],outf_sgc,logger=logger)



def ran_col_assign(randoms,data,sample_columns,tracer):
    data.rename_column('TARGETID', 'TARGETID_DATA')
    def _resamp(selregr,selregd):
        for col in sample_columns:
            randoms[col] =  np.zeros_like(data[col],shape=len(randoms))
        rand_sel = [selregr,~selregr]
        dat_sel = [ selregd,~selregd]
        for dsel,rsel in zip(dat_sel,rand_sel):
            inds = np.random.choice(len(data[dsel]),len(randoms[rsel]))
            #logger.info(str(len(data[dsel]),len(inds),np.max(inds))
            dshuf = data[dsel][inds]
            for col in sample_columns:
                randoms[col][rsel] = dshuf[col]
        
        rdl = []
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
            rdl.append(rd)
        rdr = rdl[0]/rdl[1]
        logger.info('norm factor is '+str(rdr))
        randoms['WEIGHT'][rand_sel[1]] *= rdr

    des_resamp = False
    if 'QSO' in tracer:
        des_resamp = True
    selregr = randoms['PHOTSYS'] ==  'N'
    selregd = data['PHOTSYS'] ==  'N'
    logger.info('number of randoms in N:'+str(len(randoms[selregr])))
    logger.info('number of data in N:'+str(len(data[selregd])))
    logger.info('number of randoms in S:'+str(len(randoms[~selregr])))
    logger.info('number of data in S:'+str(len(data[~selregd])))

    _resamp(selregr,selregd)
    rand_sel = [selregr,~selregr]
    dat_sel = [ selregd,~selregd]
    
    for dsel,rsel in zip(dat_sel,rand_sel):
        rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
        logger.info('data/random weighted ratio after resampling:'+str(rd))


    if des_resamp:
        logger.info('resampling in DES region')
        from regressis import footprint
        import healpy as hp
        foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_ran,phi_ran = (-randoms['DEC']+90.)*np.pi/180.,randoms['RA']*np.pi/180.
        th_dat,phi_dat = (-data['DEC']+90.)*np.pi/180.,data['RA']*np.pi/180.
        pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
        selregr = des[pixr]
        pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
        selregd = des[pixd]
        _resamp(selregr,selregd)
        rand_sel = [selregr,~selregr]
        dat_sel = [ selregd,~selregd]
    
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
            logger.info('data/random weighted ratio after resampling:'+str(rd))

    return randoms


nproc = 18


if args.tracer == 'LRG':
    zmin = 0.4
    zmax = 1.1

elif (argstracer == 'ELG_LOP') or (args.tracer == 'ELG'):
    zmin = 0.8
    zmax = 1.6

elif args.tracer == 'QSO':
    zmin = 0.8
    zmax = 2.1
elif args.tracer == 'BGS_BRIGHT-21.5':
    zmin = 0.1
    zmax = 0.4

if args.mkdat == 'y':



    selz = mock_data['Z'] > zmin
    selz &= mock_data['Z'] < zmax
    mock_data = mock_data[selz]
    logger.info('length after cutting to redshift range:'+str(len(mock_data)))
    mock_data['WEIGHT_SYS'] = np.ones(len(mock_data))
    mock_data['WEIGHT_COMP'] = np.ones(len(mock_data))
    mock_data['WEIGHT_ZFAIL'] = np.ones(len(mock_data))
    '''
    place to add imaging systematic weights and redshift failure weights would be here
    '''
    mock_data['WEIGHT'] = mock_data['WEIGHT_SYS']*mock_data['WEIGHT_COMP']*mock_data['WEIGHT_ZFAIL']
    common.write_LSS_scratchcp(mock_data,out_data_fn,logger=logger)

    #splitGC(out_data_froot,'.dat')

ran_samp_cols = ['Z','WEIGHT','WEIGHT_COMP','WEIGHT_SYS','WEIGHT_ZFAIL','TARGETID_DATA']

nran = rx-rm

randir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
ran_fname_base = randir.replace('global','dvs_ro') +'randoms-allsky-1-'

if args.mkran == 'y':
    if args.mkdat == 'n':
        mock_data = Table(fitsio.read(out_data_fn))
    def _mkran(rann):
        
        
        in_ran_fn = ran_fname_base+str(rann)+'.fits' 
        out_ran_fn = out_data_froot+str(rann)+'_clustering.ran.fits'
        rcols = ['RA','DEC']#,'PHOTSYS','TARGETID']
        common.printlog('reading random '+str(rann),logger)
        ranin = Table(fitsio.read(in_ran_fn,columns=rcols))
        common.printlog('cutting randoms '+str(rann)+' to tile area',logger)
        selY1 = is_point_in_desi(tiletab,ranin['RA'],ranin['DEC'])
        ran = ranin[selY1]
        del ranin
        logger.info(str(len(ran))+' in tiles area')
        ran = common.addNS(ran)
        ran = ran_col_assign(ran,mock_data,ran_samp_cols,args.tracer)
        common.write_LSS_scratchcp(ran,out_ran_fn,logger=logger)
        del ran
        return True
        #splitGC(out_data_froot,'.ran',rann)

    inds = np.arange(nran)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            res = pool.map(_mkran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _mkran(rn)



if args.tracer == 'QSO':
    dz = 0.02
    P0 = 6000

else:    
    dz = 0.01

if args.tracer == 'LRG':
    P0 = 10000
if args.tracer[:3] == 'ELG':
    P0 = 4000
if args.tracer[:3] == 'BGS':
    P0 = 7000


regions = ['NGC', 'SGC']

if args.nz == 'y':
    #this calculates the n(z) and then adds nbar(completeness) and FKP weights to the catalogs
    #for reg in allreg:
    fb = out_data_froot[:-1]
    fcr = fb+'_0_clustering.ran.fits'
    fcd = fb+'_clustering.dat.fits'
    fout = fb+'_nz.txt'
    common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,compmd='')
    common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,nran=nran,compmd='',par=args.par,nproc=nproc)

if args.splitGC == 'y':
    splitGC(out_data_froot,'.dat')
    def _spran(rann):
        splitGC(out_data_froot,'.ran',rann)
    inds = np.arange(nran)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            res = pool.map(_spran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _spran(rn)




