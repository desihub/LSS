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

#import LSS.mkCat_singletile.fa4lsscat as fa
#from LSS.globals import main

if os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
#parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--realization",type=int)
parser.add_argument("--prog", default="DARK")
#parser.add_argument("--veto",default='_imaging')
#parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/')
parser.add_argument("--data_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/')
parser.add_argument("--specdata_dir",help="where to find the spec data ",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is all 18)",default=18,type=int) 
parser.add_argument("--mockver", default='AbacusSummit', help = "which mocks to use. use abacus2ffa for Abacus 2nd gen fast fiber assignment")
parser.add_argument("--mockcatver", default=None, help = "if not None, gets added to the output path")

parser.add_argument("--tracer", default = 'all')
parser.add_argument("--outloc", default = None)
parser.add_argument("--par", default = 'y',help='whether to run random steps in parallel or not')
parser.add_argument("--mkdat", default = 'y')
parser.add_argument("--mkran", default = 'y')
parser.add_argument("--nz", default = 'y')
parser.add_argument("--splitGC", default = 'y')
#parser.add_argument("--remove_unassigned", default = 'y', help = 'set to n if dont want to include unassigned targets in catalog')


args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''

nside = 256
lssmapdirout = args.data_dir+'/hpmaps/'
mapn = fitsio.read(lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
maps = fitsio.read(lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
mainp = main('LRG','iron','Y1')
mapcuts = mainp.mapcuts


if args.prog == 'DARK':
    #bit = targetmask.desi_mask[args.tracer]
    bittest = targetmask.desi_mask
    desitarg='DESI_TARGET'
    if args.tracer == 'all':
        tracers = ['QSO','ELG_LOP','LRG']
    else:
        tracers = [args.tracer]
    #if args.mockver == 'abacus2ffa':
    #    tracers = [args.tracer]

#ndattot = len(mock_data)
print(tracers)


def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    #fn = Table(fitsio.read(flroot+app))
    #c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    #gc = c.transform_to('galactic')
    sel_ngc = common.splitGC(fn)#gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS(fn[sel_ngc],outf_ngc)
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS(fn[~sel_ngc],outf_sgc)



def ran_col_assign(randoms,data,sample_columns,tracer):

    def _resamp(selregr,selregd):
        for col in sample_columns:
            randoms[col] =  np.zeros(len(randoms))
        rand_sel = [selregr,~selregr]
        dat_sel = [ selregd,~selregd]
        for dsel,rsel in zip(dat_sel,rand_sel):
            inds = np.random.choice(len(data[dsel]),len(randoms[rsel]))
            print(len(data[dsel]),len(inds),np.max(inds))
            dshuf = data[dsel][inds]
            for col in sample_columns:
                randoms[col][rsel] = dshuf[col]
        randoms['WEIGHT'] *= randoms['FRAC_TLOBS_TILES'] 
        rdl = []
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
            rdl.append(rd)
        rdr = rdl[0]/rdl[1]
        print('norm factor is '+str(rdr))
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
        print('data/random weighted ratio after resampling:'+str(rd))


    if des_resamp:
        print('resampling in DES region')
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
            print('data/random weighted ratio after resampling:'+str(rd))

    return randoms

def apply_imaging_veto(ff,reccircmasks,ebits):
    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = common.maskcircandrec(ff,maskfn)
            ff = ff[~mask]

    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = common.cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))
    return ff

nproc = 9

mockdir = args.base_dir+args.mockver+'/mock'+str(args.realization)+'/'
if args.outloc == None:
    outdir = os.getenv(scratch)+'/'+args.mockver+'/mock'+str(args.realization)+'/'

if args.outloc == 'prod':
    outdir = mockdir

if args.mockcatver is not None:
    outdir += args.mockcatver + '/'

if not os.path.exists(outdir):
	os.makedirs(outdir)


print('input directory is '+mockdir)
print('output directory is '+outdir)    
    
for tracer in tracers:
   
    
    in_data_fn = mockdir + 'ffa_full_'+tracer+'.fits'

    out_data_fn = outdir+tracer+'_ffa_clustering.dat.fits'
    out_data_froot = outdir+tracer+'_ffa_'
    
    subfrac = 1
    if tracer == 'LRG':
        zmin = 0.4
        zmax = 1.1
        subfrac = 0.958 #fudge factor to get number density correct

    elif (tracer == 'ELG_LOP') or (tracer == 'ELG'):
        zmin = 0.8
        zmax = 1.6
        subfrac = .785 #determined from ration of clustering catalogs; SGC 0.77 NGC 0.793

    elif tracer == 'QSO':
        zmin = 0.8
        zmax = 2.1
        subfrac = 0.62 #determined from ratio of data with 0.8 < z < 2.1 to mock using subfrac = 1

    mainp = main(tracer,'iron','Y1')
    if args.mkdat == 'y':
    
        mock_data_tr = Table(fitsio.read(in_data_fn))
        mock_data_tr = unique(mock_data_tr,keys=['TARGETID'])
        print('length after cutting to unique targetid',len(mock_data_tr))
        selobs = mock_data_tr['WEIGHT_IIP'] != 1e20
        mock_data_tr = mock_data_tr[selobs]
        print('length after cutting to "observed" targets',len(mock_data_tr))
        mock_data_tr.rename_column('RSDZ', 'Z')
        mock_data_tr['WEIGHT_COMP'] = mock_data_tr['WEIGHT_IIP']
        
        #apply imaging vetos
        if tracer == 'LRG':
            lrgmask = fitsio.read(args.base_dir+args.mockver+'/forFA'+str(args.realization)+'_matched_input_full_lrg_imask.fits')
            mock_data_tr = join(mock_data_tr,lrgmask,keys=['TARGETID'])
            print(len(mock_data_tr))
        ebits = mainp.ebits
        reccircmasks = mainp.reccircmasks
        mock_data_tr = apply_imaging_veto(mock_data_tr,reccircmasks,ebits)
        
        mock_data_tr = common.apply_map_veto_arrays(mock_data_tr,mapn,maps,mapcuts)
        print('map data veto done')


        selz = mock_data_tr['Z'] > zmin
        selz &= mock_data_tr['Z'] < zmax
        mock_data_tr = mock_data_tr[selz]
        print('length after cutting to redshift range',len(mock_data_tr))
        sub_array = np.random.random(len(mock_data_tr))
        keep = sub_array < subfrac
        mock_data_tr = mock_data_tr[keep]
        print('length after random sub-sampling',len(mock_data_tr))
        mock_data_tr['WEIGHT_SYS'] = np.ones(len(mock_data_tr))
        mock_data_tr['WEIGHT_ZFAIL'] = np.ones(len(mock_data_tr))
        '''
        place to add imaging systematic weights and redshift failure weights would be here
        '''
        mock_data_tr['WEIGHT'] = mock_data_tr['WEIGHT_SYS']*mock_data_tr['WEIGHT_COMP']*mock_data_tr['WEIGHT_ZFAIL']
        common.write_LSS(mock_data_tr,out_data_fn)

        #splitGC(out_data_froot,'.dat')

    ran_samp_cols = ['Z','WEIGHT','WEIGHT_COMP','WEIGHT_SYS','WEIGHT_ZFAIL']

    nran = rx-rm
    ran_fname_base = args.base_dir+tracer+'_ffa_imaging_HPmapcut'

    if args.mkran == 'y':
        if args.mkdat == 'n':
            mock_data_tr = fitsio.read(out_data_fn)
        def _mkran(rann):
            
            tracerr = tracer
            if tracer == 'ELG_LOP':
                tracerr += 'notqso'
            if tracer == 'ELG':
                tracerr = 'ELG_LOPnotqso'
            in_ran_fn = ran_fname_base+str(rann)+'_full.ran.fits' 
            out_ran_fn = out_data_froot+str(rann)+'_clustering.ran.fits'
            rcols = ['RA','DEC','TILELOCID','PHOTSYS','TARGETID','NTILE','FRAC_TLOBS_TILES']
            ran = Table(fitsio.read(in_ran_fn,columns=rcols))

            ran = ran_col_assign(ran,mock_data_tr,ran_samp_cols,tracer)
            common.write_LSS(ran,out_ran_fn)
            del ran
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
        #zmin = 0.6
        #zmax = 4.5
        dz = 0.02
        P0 = 6000

    else:    
        dz = 0.01
        #zmin = 0.01
        #zmax = 1.61

    if tracer == 'LRG':
        P0 = 10000
    if tracer[:3] == 'ELG':
        P0 = 4000
    if tracer == 'BGS':
        P0 = 7000

    
    regions = ['NGC', 'SGC']

#     if args.resamp == 'y':
#         
#         for reg in regions:
#             flin = out_data_froot +reg    
#             def _parfun(rannum):
#                 ct.clusran_resamp(flin,rannum,rcols=ran_samp_cols,compmd='')#,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)
# 
#             inds = np.arange(nran)
#             if args.par == 'y':
#                 from multiprocessing import Pool
#                 with Pool(processes=nproc) as pool:
#                     res = pool.map(_parfun, inds)
#             else:
#                 for rn in range(rm,rx):
#                     _parfun(rn)

    #allreg = ['NGC', 'SGC']#'N','S',
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

    


