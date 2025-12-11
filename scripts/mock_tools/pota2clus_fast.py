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
    logger.info('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
#parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--realization",default=None)
parser.add_argument("--prog", default="DARK")
parser.add_argument("--pota_fn",help="if not None, full path to potential assignment files", default=None)
#parser.add_argument("--veto",default='_imaging')
#parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/')
parser.add_argument("--data_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/LSScats/v1/')
parser.add_argument("--in_tarf",help="full path of input target file with LRG mask info added (if it is not to be automatically found)",default=None)
#parser.add_argument("--specdata_dir",help="where to find the spec data ",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/')
parser.add_argument("--specrel",help='version of redshift catalog',default='kibo-v1')
parser.add_argument("--survey",help='the set of tiles to point to',default='DA2')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is all 18)",default=18,type=int) 
parser.add_argument("--mockver", default='AbacusSummit_v4_1', help = "which mocks to use")
parser.add_argument("--mockcatver", default=None, help = "if not None, gets added to the output path")

parser.add_argument("--tracer", default = 'all')
parser.add_argument("--outloc", default = None)
parser.add_argument("--outmd", help='write out in h5 or fits',choices=['.h5','.fits'],default = '.h5')
parser.add_argument("--par", default = 'y',help='whether to run random steps in parallel or not')
parser.add_argument("--mkdat", default = 'y')
parser.add_argument("--mkran", default = 'y')
parser.add_argument("--nz", default = 'y')
parser.add_argument("--splitGC", default = 'y')

parser.add_argument("--mk_inputran", help="this should only need to be done once and not be over-written; this script will not allow over-writes",default = 'n')
#parser.add_argument("--remove_unassigned", default = 'y', help = 'set to n if dont want to include unassigned targets in catalog')


args = parser.parse_args()
logger.info(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''

nside = 256
lssmapdirout = args.data_dir+'/hpmaps/'
if args.tracer[:3] == 'BGS':
    mapn = fitsio.read(lssmapdirout.replace('global','dvs_ro') +'BGS_BRIGHT_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout.replace('global','dvs_ro') +'BGS_BRIGHT_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')

else:
    mapn = fitsio.read(lssmapdirout.replace('global','dvs_ro') +'QSO_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout.replace('global','dvs_ro') +'QSO_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
#mainp = main('LRG','iron','Y1')



if args.tracer == 'all':
    tracers = ['QSO','LRG','ELG_LOP']
else:
    tracers = [args.tracer]

logger.info(tracers)


def read_file(fn,columns=None):
    if '.fits' in fn:
        data = Table(fitsio.read(fn.replace('global','dvs_ro')))
        if columns is not None:
            data.keep_columns(columns)
    if '.h5' in fn:
        data = common.read_hdf5_blosc(fn.replace('global','dvs_ro'),columns=columns)
    return data


def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+args.outmd
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+outmd

    fn = Table(fitsio.read(flroot.replace('global','dvs_ro') +app))
    #if datran == '.ran':
    #    fn.keep_columns(['RA', 'DEC', 'Z', 'WEIGHT', 'WEIGHT_FKP', 'TARGETID_DATA','TARGETID','NTILE'])
    #c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    #gc = c.transform_to('galactic')
    sel_ngc = common.splitGC(fn)#gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    
    outf_sgc = flroot+'SGC_'+app
    if args.outmd == '.fits':
        common.write_LSS_scratchcp(fn[sel_ngc],outf_ngc,logger=logger)
        common.write_LSS_scratchcp(fn[~sel_ngc],outf_sgc,logger=logger)

    if args.outmd == '.h5':
        common.write_LSShdf5_scratchcp(fn[sel_ngc],outf_ngc,logger=logger)
        common.write_LSShdf5_scratchcp(fn[~sel_ngc],outf_sgc,logger=logger)


def ran_col_assign(randoms,data,sample_columns,tracer,seed=0):
    
    rng = np.random.default_rng(seed=seed)
    def _resamp(selregr,selregd):
        for col in sample_columns:
            randoms[col] =  np.zeros_like(data[col],shape=len(randoms))
        rand_sel = [selregr,~selregr]
        dat_sel = [ selregd,~selregd]
        for dsel,rsel in zip(dat_sel,rand_sel):
            inds = rng.choice(len(data[dsel]),len(randoms[rsel]))
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

def apply_imaging_veto(ff,reccircmasks,ebits):
    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = common.maskcircandrec(ff,maskfn,logger=logger)
            ff = ff[~mask]

    if ebits is not None:
        logger.info('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = common.cutphotmask(ff,ebits,logger=logger)
        logger.info('number after imaging mask '+str(len(ff)))
    return ff

nproc = 18

mockdir = args.base_dir+args.mockver
if args.realization is not None:
    mockdir += '/mock'+str(args.realization)+'/'
if args.outloc == None:
    outdir = os.getenv(scratch)+'/'+args.mockver+'/mock'+str(args.realization)+'/'
elif args.outloc == 'prod':
    outdir = mockdir
else:
    outdir = args.outloc


if args.mockcatver is not None:
    outdir += args.mockcatver + '/'

if not os.path.exists(outdir):
    os.makedirs(outdir)


logger.info('input directory is '+mockdir)
logger.info('output directory is '+outdir)    

if args.pota_fn is not None:
    in_data_fn = args.pota_fn
else:
    in_data_fn = mockdir+'pota-'+args.prog+'.fits'
logger.info('using '+in_data_fn)
cols = ['LOCATION',
'FIBER',
'TARGETID',
'RA',
'DEC','RSDZ',
'PRIORITY_INIT',
'PRIORITY',
'DESI_TARGET','BRICKID','NOBS_G',
'NOBS_R',
'NOBS_Z',
'MASKBITS','ZWARN',
'COLLISION',
'TILEID']
if args.prog == 'BRIGHT':
    cols.append('R_MAG_ABS')
    mainp = main('BGS_BRIGHT',args.specrel,args.survey)

if args.prog == 'DARK':
    bittest = targetmask.desi_mask
    desitarg='DESI_TARGET'
    mainp = main('LRG',args.specrel,args.survey) #needed for bad fiber list

mapcuts = mainp.mapcuts
tsnrcut = mainp.tsnrcut
tnsrcol = mainp.tsnrcol        

if args.mkdat == 'y':
    logger.info('reading '+in_data_fn)
    mock_data = fitsio.read(in_data_fn.replace('global','dvs_ro'))#,columns=cols)
    logger.info('read '+in_data_fn.replace('global','dvs_ro'))
    selcoll = mock_data['COLLISION'] == False
    mock_data = mock_data[selcoll]
    ndattot = len(mock_data)
    lmockdat_noveto = len(mock_data)



if args.mkdat == 'y' or args.mkran == 'y':
    tilelocid = 10000*mock_data['TILEID']+mock_data['LOCATION']
    specfo =  '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.specrel+'/datcomb_'+args.prog.lower()+'_spec_zdone.fits'
    logger.info('loading specf file '+specfo)
    specf = Table(fitsio.read(specfo))
    logger.info(len(np.unique(specf['TILEID'])))
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    logger.info('loaded specf file '+specfo)
    if args.specrel == 'loa-v1' and 'v2' in args.data_dir:
        specfc = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True,logger=logger)
    else:
        specfc = common.cut_specdat(specf,badfib=mainp.badfib,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status)#common.cut_specdat(specf,badfib=mainp.badfib)
    gtl = np.unique(specfc['TILELOCID'])
    goodtl = np.isin(tilelocid,gtl)

if args.mkdat == 'y':
    mock_data = mock_data[goodtl]
    logger.info(str(lmockdat_noveto)+','+str(len(mock_data)))

    mock_data = Table(mock_data)
    mock_data = unique(mock_data,keys=['TARGETID'])
    mock_data.rename_column('RSDZ', 'Z')



    
for tracer in tracers:
    mainp = main(tracer,args.specrel,args.survey)
   
    tracerd = tracer
    #if tracer == 'BGS_BRIGHT-21.5':
    #    tracerd = 'BGS'

    out_data_fn = outdir+tracerd+'_complete_clustering.dat.'+args.outmd
    out_data_froot = outdir+tracerd+'_complete_'
    
   
    if tracer == 'LRG':
        zmin = 0.4
        zmax = 1.1

    elif (tracer == 'ELG_LOP') or (tracer == 'ELG'):
        zmin = 0.8
        zmax = 1.6

    elif tracer == 'QSO':
        zmin = 0.8
        zmax = 2.1
    elif tracer == 'BGS_BRIGHT-21.5':
        zmin = 0.1
        zmax = 0.4
    ebits = mainp.ebits
    if args.mkdat == 'y':
        if args.prog == 'DARK':
        
            bit = bittest[tracer]#targetmask.desi_mask[tracer]
            seltar = mock_data[desitarg] & bit > 0
            mock_data_tr = mock_data[seltar]
            lmockdat_noveto = len(mock_data_tr)
            logger.info('length before/after cut to target type '+tracer+' using bit '+str(bit)+' and column '+desitarg)
            logger.info(str(ndattot)+','+str(len(mock_data_tr)))
        else:
            mock_data_tr = mock_data
        lwcontam = len(mock_data_tr)
        
        sel_contam = mock_data_tr['TARGETID'] > 1e13
        mock_data_tr = mock_data_tr[~sel_contam]
        logger.info(tracer+' length before removing contaminants is '+str(lwcontam)+' and after is '+str(len(mock_data_tr)))
        if tracer == 'BGS_BRIGHT-21.5':
            if args.mockver == 'AbacusSummitBGS_v2':
                selm = (mock_data_tr['R_MAG_ABS']) < -21.5
            else:
                selm = (mock_data_tr['R_MAG_ABS']+0.05) < -21.5
            mock_data_tr = mock_data_tr[selm]
            logger.info('length after abs mag cut '+str(len(mock_data_tr)))
        
        #apply imaging vetos
        if tracer == 'LRG':
            if args.in_tarf is None:
                lrgmask = fitsio.read(args.base_dir.replace('global','dvs_ro')+args.mockver+'/forFA'+str(args.realization)+'_matched_input_full_lrg_imask.fits')
            else:
                lrgmask = fitsio.read(args.in_tarf.replace('global','dvs_ro'))
            mock_data_tr = join(mock_data_tr,lrgmask,keys=['TARGETID'])
            logger.info(len(mock_data_tr))
        
        reccircmasks = mainp.reccircmasks
        mock_data_tr = apply_imaging_veto(mock_data_tr,reccircmasks,ebits)
        
        mock_data_tr = common.apply_map_veto_arrays(mock_data_tr,mapn,maps,mapcuts)
        logger.info('map data veto done')


        selz = mock_data_tr['Z'] > zmin
        selz &= mock_data_tr['Z'] < zmax
        mock_data_tr = mock_data_tr[selz]
        logger.info('length after cutting to redshift range:'+str(len(mock_data_tr)))
        mock_data_tr['WEIGHT_SYS'] = np.ones(len(mock_data_tr))
        mock_data_tr['WEIGHT_COMP'] = np.ones(len(mock_data_tr))
        mock_data_tr['WEIGHT_ZFAIL'] = np.ones(len(mock_data_tr))
        '''
        place to add imaging systematic weights and redshift failure weights would be here
        '''
        mock_data_tr['WEIGHT'] = mock_data_tr['WEIGHT_SYS']*mock_data_tr['WEIGHT_COMP']*mock_data_tr['WEIGHT_ZFAIL']
        if args.outmd == '.fits':
            common.write_LSS_scratchcp(mock_data_tr,out_data_fn,logger=logger)
        if args.outmd == '.h5':
            common.write_LSShdf5_scratchcp(mock_data_tr,out_data_fn,logger=logger)

        #splitGC(out_data_froot,'.dat')

    ran_samp_cols = ['Z','WEIGHT','WEIGHT_COMP','WEIGHT_SYS','WEIGHT_ZFAIL','TARGETID_DATA']

    nran = rx-rm
    tracerr = tracer
    if tracer[:3] == 'BGS':
        tracerr = 'BGS_BRIGHT'
    #ran_fname_base = args.base_dir.replace('global','dvs_ro') +tracerr+'_ffa_imaging_HPmapcut'
    ran_fname_base = args.data_dir.replace('global','dvs_ro') +tracerr+'_ffa_imaging_HPmapcut'
    
    if args.mk_inputran == 'y':
        def _mk_inputran(rann):
            outfn = ran_fname_base.replace('dvs_ro','global')+str(rann)+'_full.ran.fits'
            if args.survey == 'Y1':
                infn = args.data_dir+tracerr+'_'+str(rann)+'_full_noveto.ran.fits'
            else:
                infn = args.data_dir+args.prog.lower()+'_'+str(rann)+'_full_noveto.ran.fits'
            maxp = 10000 #we don't want to apply any priority cut
        
            masked_dat = common.apply_veto(infn,outfn,ebits=ebits,zmask=False,maxp=maxp,logger=logger,reccircmasks=mainp.reccircmasks,wo='n')
            masked_dat = common.apply_map_veto_arrays(masked_dat,mapn,maps,mapcuts)
            common.write_LSS_scratchcp(masked_dat,outfn,logger=logger)
        inds = np.arange(nran)
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool(processes=9) as pool:
                res = pool.map(_mk_inputran, inds)
        else:
            for rn in inds:#range(rm,rx):
                 _mk_inputran(rn)
 

            
    if args.mkran == 'y':
        if args.mkdat == 'n':
            mock_data_tr = read_file(out_data_fn)#Table(fitsio.read(out_data_fn))
        mock_data_tr.rename_column('TARGETID', 'TARGETID_DATA')
        def _mkran(rann):
            
            tracerr = tracer
            if tracer == 'ELG_LOP':
                tracerr += 'notqso'
            if tracer == 'ELG':
                tracerr = 'ELG_LOPnotqso'
            if tracer == 'BGS_BRIGHT-21.5':
                tracerr = 'BGS_BRIGHT'
            in_ran_fn = ran_fname_base+str(rann)+'_full.ran.fits' 
            out_ran_fn = out_data_froot+str(rann)+'_clustering.ran'+args.outmd
            rcols = ['RA','DEC','PHOTSYS','TARGETID','NTILE']
            ran = Table(fitsio.read(in_ran_fn,columns=rcols))

            ran = ran_col_assign(ran,mock_data_tr,ran_samp_cols,tracer,seed=rann)
            if args.outmd == '.fits':
                common.write_LSS_scratchcp(ran,out_ran_fn,logger=logger)
            if args.outmd == '.h5':
                common.write_LSShdf5_scratchcp(ran,out_ran_fn,logger=logger)    
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
    
    

    if tracer == 'QSO':
        dz = 0.02
        P0 = 6000

    else:    
        dz = 0.01

    if tracer == 'LRG':
        P0 = 10000
    if tracer[:3] == 'ELG':
        P0 = 4000
    if tracer[:3] == 'BGS':
        P0 = 7000

    
    regions = ['NGC', 'SGC']

    if args.nz == 'y':
        #this calculates the n(z) and then adds nbar(completeness) and FKP weights to the catalogs
        #for reg in allreg:
        fb = out_data_froot[:-1]
        fcr = fb+'_0_clustering.ran'+args.outmd
        fcd = fb+'_clustering.dat'+args.outmd
        fout = fb+'_nz.txt'
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,compmd='')
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,nran=nran,compmd='',par=args.par,nproc=nproc,logger=logger,exttp=args.outmd)
    
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

    


