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
#import LSS.mkCat_singletile.fa4lsscat as fa
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
parser.add_argument("--mockver", help="type of mock to use",default='ab_firstgen')

parser.add_argument("--mockmin", help="number for the realization",default=1,type=int)
parser.add_argument("--mockmax", help="number for the realization",default=2,type=int)
parser.add_argument("--base_output", help="base directory for output",default='/pscratch/sd/a/acarnero/SecondGen/')
parser.add_argument("--targDir", help="base directory for target file",default=None)
parser.add_argument("--simName", help="base directory of AltMTL mock",default='/pscratch/sd/a/acarnero/SecondGen/altmtl_main_rea{MOCKNUM}')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--specdata", help="mountain range for spec prod",default='himalayas')
parser.add_argument("--combd", help="combine the data tiles together",default='n')
parser.add_argument("--fulld", help="make the 'full' data files ",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--add_gtl", help="whether to get the list of good tileloc from observed data",default='y')
parser.add_argument("--mkHPmaps", help="make healpix maps for imaging properties using sample randoms",default='n')
parser.add_argument("--add_veto", help="add veto column to the full files",default='n')
parser.add_argument("--apply_veto", help="apply vetos to the full files",default='n')
parser.add_argument("--apply_veto_ran", help="apply vetos to the full files",default='n')
parser.add_argument("--mkclusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--apply_map_veto", help="apply vetos to data and randoms based on values in healpix maps",default='n')
parser.add_argument("--mkclusran_allpot", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_allpot", help="make the data clustering files; these are cut to a small subset of columns",default='n')

parser.add_argument("--mkclusran_tiles", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_tiles", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--FKPfull", help="add FKP weights to full catalogs",default='n')
parser.add_argument("--split_GC",help='whether to combine N/S and then split NGC/SGC',default='n')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=1,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=2,type=int) 
parser.add_argument("--par", help="run different random number in parallel?",default='y')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--equal_data_dens", help="if y, make mock n(z) equal data n(z)", default = 'n')
parser.add_argument("--nran_clus_data", help="number of random catalogues to use for clustering data", default = 4)
parser.add_argument("--use_map_veto", help="Tag for extraveto added in name, for example, _HPmapcut", default = '')
parser.add_argument("--resamp",help="resample radial info for different selection function regions",default='n')
parser.add_argument("--getFKP", help="calculate n(z) and FKP weights on final clustering catalogs", default='n')

#--use_map_veto _HPmapcut

args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

tracer = args.tracer
survey = args.survey

if tracer[:3] == 'BGS' or tracer == 'bright' or tracer == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

randens = 10460. #the number density of randoms in the 1st gen file getting used
if args.mockver == 'ab_firstgen':
    mockdir = 'FirstGenMocks/AbacusSummit/'
    mockz = 'RSDZ'
    maindir = args.base_output +mockdir+args.survey+'/'

if args.mockver == 'EZ_3gpc1year':
    mockdir = 'FA_EZ_1year/fiberassign_EZ_3gpc/'    
    mockz = 'TRUEZ'
    maindir = args.base_output +mockdir+args.survey+'/'

if args.mockver == 'ab_secondgen':
    maindir = args.base_output
    mockz = 'RSDZ'

if args.targDir == None:
    args.targDir = maindir

if args.survey == 'MVMY1':
    tile_fn = '/global/cfs/cdirs/desi/users/FA_EZ_1year/fiberassign_EZ_3gpc/fba001/inputs/tiles.fits'
else:
    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+pr+'.fits'


tiles = fitsio.read(tile_fn)



gtl = None
if args.add_gtl == 'y':
    print('--- START ADD_GTL ---')
    datarel = args.specdata
    if args.survey == 'DA02':
        datarel = 'guadalupe'
    datadir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+datarel+'/'   
    specdat = ct.get_specdat(datadir,pdir,datarel,badfib= main(args.tracer, args.specdata, survey=args.survey).badfib)
    tlocid = 10000*specdat['TILEID'] +specdat['LOCATION']
    gtl = np.unique(tlocid)#np.unique(specdat['TILELOCID'])
    del specdat
    print('*** DONE WITH ADD_GTL ***')

def docat(mocknum, rannum):

    lssdir = os.path.join(maindir, 'mock'+str(mocknum)).format(MOCKNUM=mocknum)
    if not os.path.exists(lssdir):
        os.mkdir(lssdir)
        print('made '+lssdir)

    dirout = os.path.join(lssdir, 'LSScats')
    if not os.path.exists(dirout):
        os.mkdir(dirout)
        print('made '+dirout)


    if args.tracer != 'dark' and args.tracer != 'bright':
        if args.tracer == 'BGS_BRIGHT':
            bit = targetmask.bgs_mask[args.tracer]
            desitarg='BGS_TARGET'
        else:
            bit = targetmask.desi_mask[args.tracer]
            desitarg='DESI_TARGET'




    if args.mockver == 'ab_secondgen' and args.combd == 'y':
        print('--- START COMBD ---')
        print('entering altmtl')
        tarf = os.path.join(args.targDir, 'forFA%d.fits' % mocknum)
        ##tarf = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum #os.path.join(maindir, 'forFA_Real%d.fits' % mocknum)
        fbadir = os.path.join(args.simName, 'Univ000', 'fa', 'MAIN').format(MOCKNUM = mocknum)
        #fbadir = os.path.join(args.simName, 'Univ000', 'fa', 'MAIN').format(MOCKNUM = str(mocknum).zfill(3))
        outdir = os.path.join(maindir, 'fba' + str(mocknum)).format(MOCKNUM=mocknum)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        print('entering common.combtiles_wdup_altmtl for FASSIGN')

        asn = common.combtiles_wdup_altmtl('FASSIGN', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'assignwdup.fits'), tarf, addcols=['TARGETID','RSDZ','TRUEZ','ZWARN','PRIORITY'])
        #if using alt MTL that should have ZWARN_MTL, put that in here
        asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
        print('entering common.combtiles_wdup_altmtl for FAVAIL')
        pa = common.combtiles_wdup_altmtl('FAVAIL', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'wdup.fits'), tarf, addcols=['TARGETID','RA','DEC','PRIORITY_INIT','DESI_TARGET'])
   
        pa['TILELOCID'] = 10000*pa['TILEID'] + pa['LOCATION']
        tj = join(pa, asn, keys = ['TARGETID', 'LOCATION', 'TILEID'], join_type = 'left')
        outfs = os.path.join(lssdir, 'datcomb_' + pdir + '_tarspecwdup_zdone.fits')
        tj.write(outfs, format = 'fits', overwrite = True)
        print('wrote ' + outfs)
        tc = ct.count_tiles_better('dat', pdir, specrel = '', survey = args.survey, indir = lssdir, gtl = gtl) 
        outtc =  os.path.join(lssdir, 'Alltiles_' + pdir + '_tilelocs.dat.fits')
        tc.write(outtc, format = 'fits', overwrite = True)
        print('wrote '+outtc)
        print('*** END WITH COMBD ***')

    #specver = 'mock'    
    imbits = []   
    maxp = 3400
    if args.tracer[:3] == 'ELG':
        maxp = 3000
    if args.tracer[:3] == 'LRG' or notqso == 'notqso':
        maxp = 3200
    if args.tracer[:3] == 'BGS':
        maxp = 2100

    if args.fulld == 'y':
        print('--- START FULLD ---')
        mainp = main(args.tracer, args.specdata, survey=args.survey)
#inp        imbits = mainp.imbits
        

        #Question, is this necessary? Should be this or the list below?
        maxp = 3400
        if args.tracer[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if args.tracer[:3] == 'BGS':
            maxp = 2100

        ftar = None
        dz = os.path.join(lssdir, 'datcomb_'+pdir+'_tarspecwdup_zdone.fits')
        tlf = os.path.join(lssdir, 'Alltiles_'+pdir+'_tilelocs.dat.fits')

        fcoll = os.path.join(lssdir, 'collision_'+pdir+'_mock%d.fits' % mocknum)
        
        if not os.path.isfile(fcoll):
            fin = os.path.join(args.targDir, 'mock%d' %mocknum, 'pota-' + pr + '.fits')
            #fin = os.path.join('/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit','mock%d' %mocknum, 'pota-' + pr + '.fits')
            fcoll = mocktools.create_collision_from_pota(fin, fcoll)
        else:
            print('collision file already exist', fcoll)

        ct.mkfulldat(dz, imbits, ftar, args.tracer, bit, os.path.join(dirout, args.tracer + notqso + '_full_noveto.dat.fits'), tlf, survey = args.survey, maxp = maxp, desitarg = desitarg, specver = args.specdata, notqso = notqso, gtl_all = None, mockz = mockz,  mask_coll = fcoll, badfib = mainp.badfib, min_tsnr2 = mainp.tsnrcut, mocknum = mocknum, mockassigndir = os.path.join(args.base_output, 'fba%d' % mocknum).format(MOCKNUM=mocknum))
        print('*** END WITH FULLD ***')

#    maxp = 3400
    pthresh = 3000
    zmin = 0.8
    zmax = 3.5
    P0 = 6000
    dz_step = 0.02

    if args.tracer[:3] == 'LRG':# or notqso == 'notqso':
#        maxp = 3200
        P0 = 10000
        dz_step = 0.01
        zmin = 0.4
        zmax = 1.1
    if args.tracer[:3] == 'ELG':
        P0 = 4000
        dz_step = 0.01
#        maxp = 3000
        zmin = 0.8
        zmax = 1.6
    if args.tracer[:3] == 'BGS':
        P0 = 7000
        dz_step = 0.01
#        maxp = 2100
        pthresh = 2000
        zmin = 0.1
        zmax = 0.5
#    if notqso == 'notqso':
#        maxp = 3200

    nzmd = 'mock'
        
    if args.fullr == 'y':
        print('--- START FULLR ---')
        ldata = os.path.join(maindir, 'mock%d'% mocknum, 'datcomb_' + pdir + '_tarspecwdup_zdone.fits').format(MOCKNUM=mocknum)
        specft = fitsio.read(ldata) #Is this from data or mock? 
        wg = np.isin(specft['TILELOCID'], gtl)
        specft = Table(specft[wg])
        lznp = common.find_znotposs(specft)
        del specft
        mainp = main(args.tracer, args.specdata, survey=args.survey)
        imbits = mainp.imbits
        tsnrcut = mainp.tsnrcut
        global _parfun1
        def _parfun1(rann):
            ranfile = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'rancomb_%d%swdupspec_zdone.fits' % (rann, pdir)) 
            alltileloc = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'rancomb_%d%s_Alltilelocinfo.fits' % (rann, pdir)) 
                
            ranfile, alltileloc = mocktools.createrancomb_wdupspec(lssdir, ranfile, alltileloc, os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup.fits').format(MOCKNUM=mocknum), os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'datcomb_'+pdir+'_spec_zdone.fits'))
            outf = os.path.join(dirout, args.tracer+notqso+'_'+str(rann)+'_full_noveto.ran.fits')
            ct.mkfullran(gtl, lznp, os.path.join(maindir, 'mock'+str(mocknum)).format(MOCKNUM=mocknum), rann, imbits, outf, args.tracer, pdir, notqso = notqso, maxp = maxp, min_tsnr2 = tsnrcut) 
##        ct.mkfullran(gtlf,lznp,lssdir,rannum,imbits,outf,args.tracer,pdir,notqso=notqso,maxp=maxp,tlid_full=tlid_full)
        if args.par == 'n':
            for rn in range(rannum[0], rannum[1]):
                _parfun1(rn)
        else:
            from multiprocessing import Pool

            inds = np.arange(rannum[0], rannum[1])
            (rannum[1]-rannum[0])*2
            nproc = 9 #try this so doesn't run out of memory
            with Pool(processes=nproc) as pool:
                res = pool.map(_parfun1, inds)
                pool.close()
                pool.join()
        print('*** END WITH FULLR ***')

    mainp = main(args.tracer, args.specdata, survey=args.survey)

    if args.apply_veto == 'y':
        print('--- START APPLY_VETO ---')
        print('applying vetos to mock ' + str(mocknum))
        fin = os.path.join(dirout, args.tracer + notqso + '_full_noveto.dat.fits')
        if args.tracer == 'LRG' and args.survey == 'Y1' and args.mockver == 'ab_secondgen':
            lrgf = Table.read(fin)
            if 'lrg_mask' not in lrgf.columns:
                
                lrgmask = Table.read(os.path.join(args.targDir, 'forFA%d_matched_input_full_lrg_imask.fits' % mocknum))
                lrgf = join(lrgf, lrgmask, keys=['TARGETID'])
                common.write_LSS(lrgf, fin)

        elif args.survey=='Y1' and args.mockver == 'ab_secondgen':
            novf = Table.read(fin)
            maskcols = []
            cols = ['NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']
            ovw = False
            for colm in cols:
                if colm not in novf.columns:
                    maskcols.append(colm)
            if len(maskcols) > 0:
                maskcols.append('TARGETID')
                targf = Table(fitsio.read(os.path.join(args.targDir, 'forFA%d.fits' % mocknum), columns = maskcols))
                novf = join(novf, targf, keys=['TARGETID'])
                ovw = True
            if 'PHOTSYS' not in novf.columns:
                novf['PHOTSYS'] = 'N'
                sel = novf['DEC'] < 32.375
                wra = (novf['RA'] > 100-novf['DEC'])
                wra &= (novf['RA'] < 280 +novf['DEC'])
                sel |= ~wra
                novf['PHOTSYS'][sel] = 'S'
                ovw = True
            if ovw:
                common.write_LSS(novf, fin)

        fout = os.path.join(dirout, args.tracer + notqso + '_full.dat.fits')
        common.apply_veto(fin, fout, ebits = mainp.ebits, zmask = False, maxp = maxp, reccircmasks = mainp.reccircmasks)
        
        global _parfun2
        def _parfun2(rann):
            print('applying vetos to random ' + str(rann))
            fin = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full_noveto.ran.fits')
            fout = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full.ran.fits')
            if args.tracer != 'QSO':
                common.add_veto_col(fin, ran = True, tracer_mask = args.tracer[:3].lower(), rann = rann)
            common.apply_veto(fin, fout, ebits = mainp.ebits, zmask = False, maxp = maxp, reccircmasks = mainp.reccircmasks)
        
        if args.par == 'n':
            for rn in range(rannum[0], rannum[1]):
                _parfun2(rn)
        else:
            from multiprocessing import Pool

            inds = np.arange(rannum[0], rannum[1])
            (rannum[1]-rannum[0])*2
            nproc = 9 #try this so doesn't run out of memory
            with Pool(processes=nproc) as pool:
                res = pool.map(_parfun2, inds)
                pool.close()
                pool.join()
        print('*** END VETO ***')
        #print('random veto '+str(ii)+' done')
    if args.apply_veto_ran == 'y':
        for rann in range(rannum[0], rannum[1]):
            fin = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full_noveto.ran.fits')
            fout = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full.ran.fits')
            if args.tracer != 'QSO':
                common.add_veto_col(fin, ran = True, tracer_mask = args.tracer[:3].lower(), rann = rann)
            common.apply_veto(fin, fout, ebits = mainp.ebits, zmask = False, maxp = maxp, reccircmasks = mainp.reccircmasks)

    regl = ['_N','_S']    


    wzm = ''

    lssmapdirout = os.path.join(dirout, 'hpmaps')
    nside = 256
    tracer_clus = args.tracer + notqso + wzm

    if args.mkHPmaps == 'y':
        print('--- START MKHPMAPS ---')
        from LSS.imaging.sky_maps import create_pixweight_file, rancat_names_to_pixweight_name
        if not os.path.exists(lssmapdirout):
            os.mkdir(lssmapdirout)
            print('made ' + lssmapdirout)
        new_cols=mainp.new_cols#['STARDENS','HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z', 'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15', 'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK']
        fid_cols=mainp.fid_cols#['EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
        


        fieldslist = new_cols + fid_cols
        lssmapdir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/external_input_maps/'
        rancatname = os.path.join(dirout, tracer_clus + '_*_full.ran.fits')
        rancatlist = sorted(glob.glob(rancatname))
        for temp_rannum in range(0, 18): #rancat in rancatlist:

            rancat = os.path.join(dirout, tracer_clus + '_%d_full.ran.fits' % temp_rannum) 
            
            print('filling randoms with imaging properties for', rancat, temp_rannum)
            common.add_map_cols(rancat, temp_rannum, new_cols = new_cols, fid_cols = fid_cols)
            print('done with ', str(temp_rannum))

        masklist = list(np.zeros(len(fieldslist), dtype = int))
        for reg in ['N','S']:
            outfn = os.path.join(lssmapdirout, tracer_clus + '_mapprops_healpix_nested_nside' + str(nside) + '_' + reg + '.fits')
            create_pixweight_file(rancatlist, fieldslist, masklist, nside_out = nside,
                          lssmapdir = lssmapdir, outfn = outfn, reg = reg)    
        print('*** END WITH MKHPMAPS ***')

    if args.apply_map_veto == 'y':
        print('--- START APPLY_MAP_VETO ---')
        import healpy as hp
        lssmapdirout = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/hpmaps'
        mapn = fitsio.read(os.path.join(lssmapdirout, tracer_clus + '_mapprops_healpix_nested_nside' + str(nside) + '_N.fits'))
        maps = fitsio.read(os.path.join(lssmapdirout, tracer_clus + '_mapprops_healpix_nested_nside' + str(nside) + '_S.fits'))
        mapcuts = mainp.mapcuts

#        if args.ranonly != 'y':
        fout = os.path.join(dirout, args.tracer + notqso + '_full.dat.fits')
        fin = fout.replace('global', 'dvs_ro')
        fout = fout.replace('_full', '_full_HPmapcut')
        common.apply_map_veto(fin, fout, mapn, maps, mapcuts)
        print('data veto done, now doing randoms')
        
        global _parfun3
        def _parfun3(rn):
            fout = os.path.join(dirout, args.tracer + notqso + '_' + str(rn) + '_full.ran.fits')
            fin = fout.replace('global', 'dvs_ro')
            fout = fout.replace('_full', '_full_HPmapcut')
            common.apply_map_veto(fin, fout, mapn, maps, mapcuts)
            print('random veto ' + str(rn) + ' done')
        if args.par == 'n':
            for rn in range(rannum[0], rannum[1]):
                _parfun3(rn)
        else:
            from multiprocessing import Pool

            inds = np.arange(rannum[0], rannum[1])
            (rannum[1] - rannum[0])*2
            nproc = 9 #try this so doesn't run out of memory
            with Pool(processes = nproc) as pool:
                res = pool.map(_parfun3, inds)
                pool.close()
                pool.join()
        print('*** END WITH APPLY_MAP_VETO ***')

    nztl = []
    if args.mkclusdat == 'y':
        print('--- START MKCLUSDAT ---')
        nztl.append('')
        fin = os.path.join(dirout, args.tracer + notqso + '_full' + args.use_map_veto + '.dat.fits')
        #ct.mkclusdat(os.path.join(dirout,args.tracer+notqso),tp=args.tracer,dchi2=None,tsnrcut=0,zmin=zmin,zmax=zmax)#,ntilecut=ntile)
        ct.mkclusdat(os.path.join(dirout, args.tracer + notqso), tp = args.tracer, dchi2 = None, tsnrcut = 0, zmin = zmin, zmax = zmax, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
        #ct.mkclusdat(os.path.join(dirout, args.tracer + notqso), tp = args.tracer, dchi2 = None, splitNS='y', tsnrcut = 0, zmin = zmin, zmax = zmax, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
        print('*** END WITH MKCLUSDAT ***')

    
    if args.equal_data_dens == 'y':
        data_dir = "/dvs_ro/project/projectdirs/desi/survey/catalogs/edav1/da02/LSScats/clustering"
        
        for i in ['N','S']:
            mocktools.mock_equal_data_density(dirout, data_dir, dirout, args.tracer, i, zmin, zmax, args.nran_clus_data, randens)
        
        

    if args.mkclusran == 'y':
        print('--- START MKCLUSRAN ---')
        if len(nztl) == 0:
            nztl.append('')
        rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        fl = os.path.join(dirout, args.tracer + notqso + '_')
        print('adding tlobs to randoms with ', fl)
        
        global _parfun4
        def _parfun4(rann):
            ct.add_tlobs_ran(fl, rann, hpmapcut = args.use_map_veto)
            ct.mkclusran(os.path.join(dirout, args.tracer + notqso + '_'), os.path.join(dirout, args.tracer + notqso + '_'), rann, rcols = rcols,  tsnrcut = 0, tsnrcol = tsnrcol, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
            #ct.mkclusran(os.path.join(dirout, args.tracer + notqso + '_'), os.path.join(dirout, args.tracer + notqso + '_'), rann, rcols = rcols, nosplit='n', tsnrcut = 0, tsnrcol = tsnrcol, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        if args.par == 'n':
            for rn in range(rannum[0], rannum[1]):
                _parfun4(rn)
        else:
            from multiprocessing import Pool

            inds = np.arange(rannum[0], rannum[1])
            (rannum[1]-rannum[0])*2
            nproc = 9 #try this so doesn't run out of memory
            with Pool(processes=nproc) as pool:
                res = pool.map(_parfun4, inds)
                pool.close()
                pool.join()

#        ct.clusNStoGC(os.path.join(dirout, args.tracer + notqso+'_'), rannum[1] - rannum[0])
        print('*** END WITH MKCLUSRAN ***')

    nran = rannum[1] - rannum[0]
##    regions = ['NGC', 'SGC']#, 'N', 'S']


    regions = ['NGC', 'SGC', 'N', 'S']

    if args.resamp == 'y':
        ct.splitclusGC(os.path.join(dirout, args.tracer + notqso + '_'), nran)
        ct.splitclusNS(os.path.join(dirout, args.tracer + notqso + '_'), nran)
        print('--- START RESAMP AND NZ ---')
        rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']#, 'WEIGHT_FKP']
        
        for reg in regions:
            fb = os.path.join(dirout, tracer_clus + '_' + reg)
            
            global _parfun5
            def _parfun5(rann):
                ct.clusran_resamp(fb, rann, rcols = rcols)#,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)

            inds = np.arange(rannum[0], rannum[1])
            if args.par == 'y':
                from multiprocessing import Pool

                with Pool(processes = nran*2) as pool:
                    res = pool.map(_parfun5, inds)
                    pool.close()
                    pool.join()
            else:
                for rn in range(rannum[0], rannum[1]):
                    _parfun5(rn)

            fcr = fb + '_' + str(rannum[0]) + '_clustering.ran.fits'
            fcd = fb + '_clustering.dat.fits'
            fout = fb + '_nz.txt'
            common.mknz(fcd, fcr, fout, bs = dz_step, zmin = zmin, zmax = zmax)
            common.addnbar(fb, bs = dz_step, zmin = zmin, zmax = zmax, P0 = P0, nran = nran)
        print('*** END WITH RESAMP AND NZ ***')

    '''
    if args.FKPfull == 'y':
    
        fb = os.path.join(dirout, tracer_clus)
        fbr = fb
        if type == 'BGS_BRIGHT-21.5':
            fbr = dirout+'BGS_BRIGHT'

        fcr = fbr +' _0_full'+args.use_map_veto+'.ran.fits'
        fcd = fb + '_full'+args.use_map_veto+'.dat.fits'
        nz = common.mknz_full(fcd, fcr, args.tracer[:3], bs=dz_step, zmin=zmin, zmax=zmax)
        common.addFKPfull(fcd, nz, args.tracer[:3], bs=dz, zmin=zmin, zmax=zmax, P0=P0)
    '''

    if args.mkclusdat_allpot == 'y':
        #fbadir = maindir+'fba'+str(mocknum)
        #tarf = fbadir+'/targs.fits'
        #'/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit
        tarf = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum
        #args.base_output +mockdir+'/forFA'+str(mocknum)+'.fits'
        ztab = Table(fitsio.read(tarf, columns=['TARGETID','RSDZ']))
        ztab.rename_column('RSDZ', 'Z')
        
        mocktools.mkclusdat_allpot(os.path.join(dirout, args.tracer + notqso), ztab, tp = args.tracer, dchi2 = None, tsnrcut = 0, zmin = zmin, zmax = zmax)#,ntilecut=ntile)



    if args.mkclusran_allpot == 'y':
        nztl.append('_complete')
        rcols=['Z','WEIGHT']
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        fl = os.path.join(dirout, args.tracer + notqso + '_')
        for rann in range(rannum[0], rannum[1]):
            ct.add_tlobs_ran(fl, rann)
            ct.mkclusran(os.path.join(dirout, args.tracer + notqso + '_'), os.path.join(dirout, args.tracer + notqso + '_complete' + '_'), rann, rcols = rcols, tsnrcut = 0, tsnrcol = tsnrcol, compmd = '')#,ntilecut=ntile,ccut=ccut)

        ct.clusNStoGC(os.path.join(dirout, args.tracer + notqso + '_complete_'), nran=1)
    if args.mkclusdat_tiles == 'y':
        fbadir = maindir + 'fba' + str(mocknum)
        tarf = fbadir + '/targs.fits'

        ztab = Table(fitsio.read(tarf, columns=['RA','DEC','RSDZ','DESI_TARGET']))
        ztab.rename_column('RSDZ', 'Z')
        mocktools.mkclusdat_tiles(dirout + args.tracer + notqso, ztab, bit, zmin = zmin, zmax = zmax)#,ntilecut=ntile)

    if args.mkclusran_tiles == 'y':
        nztl.append('_tiles')
        fbadir = maindir+'random_fba'+str(rannum)
        tarf = fbadir+'/targs.fits'
        ffc = Table(fitsio.read(tarf,columns=['RA','DEC']))
        rcols=['Z','WEIGHT']
        mocktools.mkclusran_tiles(ffc,dirout+args.tracer+notqso+'_tiles'+'_',rannum,rcols=rcols)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        for reg in regl:
            ranf = dirout+args.tracer+notqso+'_tiles'+reg+'_'+str(rannum)+'_clustering.ran.fits'
            ranfm = dirout+args.tracer+notqso+'_tiles'+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
            os.system('mv '+ranf+' '+ranfm)

   

    #print('done with random '+str(ii))
    return True
        #ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if __name__ == '__main__':
    rx = args.maxr
    rm = args.minr
    mockmin = args.mockmin
    mockmax = args.mockmax
    for i in range(mockmin, mockmax):
        docat(i, [rm, rx])
