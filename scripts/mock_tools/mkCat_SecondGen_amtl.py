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
from astropy.table import Table,join,unique,vstack,setdiff
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi


import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
#import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main
import errno

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made ' + value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--mockver", help="type of mock to use",default='ab_firstgen')

parser.add_argument("--mocknum", help="number for the realization",default=1,type=int)
parser.add_argument("--ccut", help="extra-cut",default=None)
parser.add_argument("--base_output", help="base directory for output",default=os.getenv('SCRATCH')+'/SecondGen/')
parser.add_argument("--outmd", help="whether to write in scratch",default='scratch')
parser.add_argument("--targDir", help="base directory for target file",default=None)
parser.add_argument("--simName", help="base directory of AltMTL mock",default='/pscratch/sd/a/acarnero/SecondGen/altmtl_main_rea{MOCKNUM}')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--specdata", help="mountain range for spec prod",default='iron')
parser.add_argument("--combd", help="combine the data tiles together",default='n')
parser.add_argument("--joindspec", help="combine the target and spec info together",default='n')
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

parser.add_argument("--start_from_full",help="whether to start from the full catalogs already moved the the final directory",default='n')

parser.add_argument("--mkclusran_tiles", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_tiles", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--FKPfull", help="add FKP weights to full catalogs",default='n')
parser.add_argument("--splitGC",help='whether to combine N/S and then split NGC/SGC',default='n')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=1,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=2,type=int) 
parser.add_argument("--par", help="run different random number in parallel?",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--equal_data_dens", help="if y, make mock n(z) equal data n(z)", default = 'n')
parser.add_argument("--nran_clus_data", help="number of random catalogues to use for clustering data", default = 4)
parser.add_argument("--use_map_veto", help="Tag for extraveto added in name, for example, _HPmapcut", default = '')
parser.add_argument("--resamp",help="resample radial info for different selection function regions",default='n')
parser.add_argument("--getFKP", help="calculate n(z) and FKP weights on final clustering catalogs", default='n')
parser.add_argument("--add_bitweights", help="Add bitweights to files before creating the final clustering catalogs.", default=None)
parser.add_argument("--add_weight_ntile", help="Add NTILE weights to full catalogs to make it compatible with PIP and angular upweithing", default='n')

#--use_map_veto _HPmapcut

import logging

# create logger
logname = 'LSSran'
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


args = parser.parse_args()
print(args)

mocknum = args.mocknum

rm = int(args.minr)
rx = int(args.maxr)
rannum = (rm,rx)

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


maindir = args.base_output
mockz = 'RSDZ'

if args.targDir == None:
    args.targDir = maindir


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


lssdir = os.path.join(maindir, 'mock'+str(mocknum)).format(MOCKNUM=mocknum)
test_dir(lssdir)
#if not os.path.exists(lssdir):
#    os.mkdir(lssdir)
#    print('made '+lssdir)

dirout = os.path.join(lssdir, 'LSScats')
dirfinal = dirout
if args.outmd == 'scratch':
    dirout = dirout.replace('/global/cfs/cdirs/desi/survey/catalogs/',os.getenv('SCRATCH')+'/')
test_dir(dirout)

#if not os.path.exists(dirout):
#    os.makedirs(dirout)
#    print('made '+dirout)


if args.tracer != 'dark' and args.tracer != 'bright':
    if args.tracer == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[args.tracer]
        #desitarg='DESI_TARGET'
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[args.tracer]
        desitarg='DESI_TARGET'



asn = None
pa = None
outdir = os.path.join(maindir, 'fba' + str(mocknum)).format(MOCKNUM=mocknum)
test_dir(outdir)

if args.mockver == 'ab_secondgen' and args.combd == 'y':
    print('--- START COMBD ---')
    print('entering altmtl')
    tarf = os.path.join(args.targDir, 'forFA%d.fits' % mocknum)
    ##tarf = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum #os.path.join(maindir, 'forFA_Real%d.fits' % mocknum)
    fbadir = os.path.join(maindir, 'Univ000', 'fa', 'MAIN').format(MOCKNUM = mocknum)
    #fbadir = os.path.join(args.simName, 'Univ000', 'fa', 'MAIN').format(MOCKNUM = str(mocknum).zfill(3))
    print('entering common.combtiles_wdup_altmtl for FASSIGN')

    asn = common.combtiles_wdup_altmtl('FASSIGN', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'assignwdup.fits'), tarf, addcols=['TARGETID','RSDZ','TRUEZ','ZWARN'])
    #asn = common.combtiles_assign_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RSDZ','TRUEZ','ZWARN'],fba=True,tp='dark')
    #if using alt MTL that should have ZWARN_MTL, put that in here
    asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
    print('entering common.combtiles_wdup_altmtl for FAVAIL')

    cols = ['TARGETID','RA','DEC','PRIORITY_INIT','DESI_TARGET']
    if pdir == 'bright':
        cols.append('BGS_TARGET')
        cols.append('R_MAG_ABS')
        cols.append('G_R_OBS')
        cols.append('G_R_REST')
    pa = common.combtiles_wdup_altmtl('FAVAIL', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'wdup.fits'), tarf, addcols=cols)

fcoll = os.path.join(lssdir, 'collision_'+pdir+'_mock%d.fits' % mocknum)
if args.joindspec == 'y':

    if asn is None:
        afn = os.path.join(outdir, 'datcomb_' + pdir + 'assignwdup.fits')
        asn = fitsio.read(afn)
        print('loaded assignments')
    if pa is None:
        pafn = os.path.join(outdir, 'datcomb_' + pdir + 'wdup.fits')
        pa = Table(fitsio.read(pafn))
        print('loaded potential assignements')
    pa['TILELOCID'] = 10000*pa['TILEID'] + pa['LOCATION']
    print('about to join assignments and potential assignments')
    tj = join(pa, asn, keys = ['TARGETID', 'LOCATION', 'TILEID'], join_type = 'left')
    
    
    if not os.path.isfile(fcoll):
        fin = os.path.join(args.targDir, 'mock%d' %mocknum, 'pota-' + pr + '.fits')
        #fin = os.path.join('/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit','mock%d' %mocknum, 'pota-' + pr + '.fits')
        fcoll = mocktools.create_collision_from_pota(fin, fcoll)
    else:
        print('collision file already exist', fcoll)

    coll = Table(fitsio.read(fcoll))
    print('length before masking collisions '+str(len(tj)))
    tj = setdiff(tj,coll,keys=['TARGETID','LOCATION','TILEID'])
    print('length after masking collisions '+str(len(tj)))

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

    ftar = None
    dz = os.path.join(lssdir, 'datcomb_'+pdir+'_tarspecwdup_zdone.fits')
    tlf = os.path.join(lssdir, 'Alltiles_'+pdir+'_tilelocs.dat.fits')

    ct.mkfulldat_mock(dz, imbits, ftar, args.tracer, bit, os.path.join(dirout, args.tracer + notqso + '_full_noveto.dat.fits'), tlf, survey = args.survey, maxp = maxp, desitarg = desitarg, specver = args.specdata, notqso = notqso, gtl_all = None, mockz = mockz,  mask_coll = fcoll, badfib = mainp.badfib, min_tsnr2 = mainp.tsnrcut, mocknum = mocknum, mockassigndir = os.path.join(args.base_output, 'fba%d' % mocknum).format(MOCKNUM=mocknum))
    print('*** END WITH FULLD ***')

#    maxp = 3400
pthresh = 3000
zmin = 0.8
zmax = 3.5
P0 = 6000
dz_step = 0.02

zsplit = None
subfrac = 1
if tracer == 'QSO':
    zmin = 0.8
    zmax = 2.1
    subfrac = 0.66 #determined from ratio of data with 0.8 < z < 2.1 to mock using subfrac = 1 for altmtl version 3_1

if args.tracer[:3] == 'LRG':# or notqso == 'notqso':
#        maxp = 3200
    P0 = 10000
    dz_step = 0.01
    zmin = 0.4
    zmax = 1.1
    subfrac = 0.976
if args.tracer[:3] == 'ELG':
    P0 = 4000
    dz_step = 0.01
#        maxp = 3000
    zmin = 0.8
    zmax = 1.6
    subfrac = [0.69,0.54]#0.676
    zsplit=1.5
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
    #lznp = common.find_znotposs(specft) #doesn't actually get used and takes a long time
    lznp = None
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
        #(rannum[1]-rannum[0])*2
        nproc = rx-rm #try 9 if runs out of memory
        ####HERE nproc = 18 #try 9 if runs out of memory
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun1, inds)
            pool.close()
            pool.join()
    print('*** END WITH FULLR ***')

mainp = main(args.tracer, args.specdata, survey=args.survey)

tracer_clus = args.tracer + notqso 

if args.apply_veto == 'y':
    print('--- START APPLY_VETO; including HP maps---')
    print('applying vetos to mock ' + str(mocknum))
    import healpy as hp
    nside = 256
    lssmapdirout = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/hpmaps'
    mapn = fitsio.read(os.path.join(lssmapdirout, tracer_clus + '_mapprops_healpix_nested_nside' + str(nside) + '_N.fits'))
    maps = fitsio.read(os.path.join(lssmapdirout, tracer_clus + '_mapprops_healpix_nested_nside' + str(nside) + '_S.fits'))
    mapcuts = mainp.mapcuts

    fin = os.path.join(dirout, args.tracer + notqso + '_full_noveto.dat.fits')
    colnames = list(fitsio.read(fin,rows=1).dtype.names)
    maskcols = ['NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']
    addlrg = 0
    if args.tracer == 'LRG':
        if 'lrg_mask' not in colnames:
            addcols = 1
            addlrg = 1
    coltest = np.isin(maskcols,colnames)
    readcols = maskcols.copy()
    readcols.append('TARGETID')
    addcols = 0
    if np.sum(coltest) != len(maskcols):
        addcols = 1
        joinmask = 1
        print(maskcols,coltest,colnames)
        #sys.exit()
    if 'PHOTSYS' not in colnames:
        addcols = 1
    if addcols == 1:
        dataf = Table(fitsio.read(fin))
        if addlrg == 1:
            lrgmask = Table.read(os.path.join(args.targDir.replace('global','dvs_ro'), 'forFA%d_matched_input_full_lrg_imask.fits' % mocknum))
            dataf = join(dataf, lrgmask, keys=['TARGETID'])
        if joinmask == 1:                   
            targf = Table(fitsio.read(os.path.join(args.targDir.replace('global','dvs_ro'), 'forFA%d.fits' % mocknum), columns = readcols))
            dataf = join(dataf, targf, keys=['TARGETID'])
        if 'PHOTSYS' not in colnames:
            dataf = common.addNS(dataf)
        common.write_LSS(dataf, fin)


    fout = os.path.join(dirout, args.tracer + notqso + '_full'+args.use_map_veto + '.dat.fits')
    dataf = common.apply_veto(fin, fout,ebits = mainp.ebits, zmask = False, maxp = maxp, reccircmasks = mainp.reccircmasks,wo='n',mapveto=args.use_map_veto) #returns vetoed array
    dataf = common.apply_map_veto_arrays(dataf,mapn,maps,mapcuts)
    common.write_LSS(dataf,fout)
    print('data veto done, now doing randoms')
    
    global _parfun2
    def _parfun2(rann):
        print('applying vetos to random ' + str(rann))
        fin = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full_noveto.ran.fits')
        fout = os.path.join(dirout, args.tracer + notqso + '_' + str(rann) + '_full'+args.use_map_veto + '.ran.fits')
        if args.tracer == 'LRG':
            common.add_veto_col(fin, ran = True, tracer_mask = args.tracer[:3].lower(), rann = rann)
        ranf = common.apply_veto(fin, ebits = mainp.ebits, zmask = False, maxp = maxp, reccircmasks = mainp.reccircmasks)
        ranf = common.apply_map_veto_arrays(ranf,mapn,maps,mapcuts)
        common.write_LSS(ranf,fout)
    
    if args.par == 'n':
        for rn in range(rannum[0], rannum[1]):
            _parfun2(rn)
    else:
        from multiprocessing import Pool

        inds = np.arange(rannum[0], rannum[1])
        nproc = 18 #try this so doesn't run out of memory
        if tracer == 'QSO' or tracer == 'LRG':
            nproc = 9 #QSO has OOM with all 18
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun2, inds)
            pool.close()
            pool.join()
    print('*** END VETO ***')
    #print('random veto '+str(ii)+' done')

readdir = dirout
if args.start_from_full == 'y':
	readdir = dirfinal


nztl = []
if args.mkclusdat == 'y':
    print('--- START MKCLUSDAT ---')
    nztl.append('')
    
    if args.ccut is not None:
        ffile = Table.read(os.path.join(readdir, args.tracer + notqso + '_full'+args.use_map_veto + '.dat.fits').replace('global','dvs_ro'))
        if 'R_MAG_ABS' not in ffile.columns:
            targets = Table(fitsio.read(os.path.join(args.targDir, 'forFA{MOCKNUM}.fits').format(MOCKNUM=mocknum).replace('global','dvs_ro'), columns=['TARGETID', 'R_MAG_ABS']))
            nm = Table(join(ffile, targets, keys=['TARGETID']))
        #print(nm)
            common.write_LSS(nm, os.path.join(readdir, args.tracer + notqso + '_full'+args.use_map_veto + '.dat.fits'))
        #nm.write(ffile, overwrite=True)

    if args.add_bitweights is not None:
        ffile = Table.read(os.path.join(readdir, args.tracer + notqso + '_full'+args.use_map_veto + '.dat.fits').replace('global','dvs_ro'))
        bitweights_file = Table.read(args.add_bitweights)
        nm = Table(join(ffile, bitweights_file, join_type='left', keys=['TARGETID']))
        common.write_LSS(nm, os.path.join(readdir, args.tracer + notqso + '_full'+args.use_map_veto + '.dat.fits'))

    ct.mkclusdat(os.path.join(readdir, args.tracer + notqso), tp = args.tracer, dchi2 = None, tsnrcut = 0, zmin = zmin, zmax = zmax, use_map_veto = args.use_map_veto, subfrac=subfrac,zsplit=zsplit, ismock=True, ccut=args.ccut)#,ntilecut=ntile,ccut=ccut)
    #ct.mkclusdat(os.path.join(dirout, args.tracer + notqso), tp = args.tracer, dchi2 = None, splitNS='y', tsnrcut = 0, zmin = zmin, zmax = zmax, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
    print('*** END WITH MKCLUSDAT ***')


    
    
finaltracer = args.tracer + notqso #+ '_'
if args.tracer[:3] == 'BGS':
    if args.ccut is not None:
        finaltracer = args.tracer + str(args.ccut) #+ '_'

rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL','TARGETID_DATA']
if args.mkclusran == 'y':
    print('--- START MKCLUSRAN ---')
    if len(nztl) == 0:
        nztl.append('')

    tsnrcol = 'TSNR2_ELG'
    if args.tracer[:3] == 'BGS':
        fl = os.path.join(dirout, finaltracer) + '_'
        cols_clustering = Table.read(fl.replace('global','dvs_ro')+'clustering.dat.fits').columns
        if 'G_R_OBS' in cols_clustering:
            rcols.append('G_R_OBS')
        if 'G_R_REST' in cols_clustering:
            rcols.append('G_R_REST')
        if 'R_MAG_ABS' in cols_clustering:
            rcols.append('R_MAG_ABS')

        tsnrcol = 'TSNR2_BGS'
        if args.ccut is not None:
            for rn in range(rannum[0], rannum[1]):
                if not os.path.isfile('%s%s_%d_full_HPmapcut.ran.fits'% (os.path.join(dirout, args.tracer), str(args.ccut), rn)):
                    os.system('cp %s_%d_full_HPmapcut.ran.fits %s%s_%d_full_HPmapcut.ran.fits' %(os.path.join(dirout, args.tracer), rn, os.path.join(dirout, args.tracer), str(args.ccut), rn))
                #print('cp %s_%d_full_HPmapcut.ran.fits %s%s_%d_full_HPmapcut.ran.fits' %(os.path.join(dirout, args.tracer), rn, os.path.join(dirout, args.tracer), str(args.ccut), rn))
            os.system('cp %s_frac_tlobs.fits %s%s_frac_tlobs.fits' %(os.path.join(dirout, args.tracer), os.path.join(dirout, args.tracer), str(args.ccut)))
    fl = os.path.join(dirout, finaltracer) + '_'
    print('adding tlobs to randoms with ', fl)
    clus_arrays = [fitsio.read(fl.replace('global','dvs_ro')+'clustering.dat.fits')]
    global _parfun4
    def _parfun4(rann):
        #ct.add_tlobs_ran(fl, rann, hpmapcut = args.use_map_veto)
#        print(os.path.join(readdir, finaltracer) + '_', os.path.join(dirout, finaltracer) + '_', rann, rcols, -1, tsnrcol, args.use_map_veto,  clus_arrays, 'y')
        ct.mkclusran(os.path.join(readdir, finaltracer) + '_', os.path.join(dirout, finaltracer) + '_', rann, rcols = rcols,  tsnrcut = -1, tsnrcol = tsnrcol, use_map_veto = args.use_map_veto,clus_arrays=clus_arrays,add_tlobs='y')#,ntilecut=ntile,ccut=ccut)
        #ct.mkclusran(os.path.join(dirout, args.tracer + notqso + '_'), os.path.join(dirout, args.tracer + notqso + '_'), rann, rcols = rcols, nosplit='n', tsnrcut = 0, tsnrcol = tsnrcol, use_map_veto = args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
    #for clustering, make rannum start from 0
    if args.par == 'n':
        for rn in range(rannum[0], rannum[1]):
            _parfun4(rn)
    else:
        from multiprocessing import Pool

        inds = np.arange(rannum[0], rannum[1])
        nproc = 18 #try this so doesn't run out of memory
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun4, inds)
            pool.close()
            pool.join()

#        ct.clusNStoGC(os.path.join(dirout, args.tracer + notqso+'_'), rannum[1] - rannum[0])
    print('*** END WITH MKCLUSRAN ***')


fb = os.path.join(dirout, finaltracer)
nran = rx-rm



nproc = 18
##fb = os.path.join(dirout, tracer_clus)
if args.nz == 'y':
    #this calculates the n(z) and then adds nbar(completeness) and FKP weights to the catalogs
    #for reg in allreg:
    
    fcr = fb+'_0_clustering.ran.fits'
    fcd = fb+'_clustering.dat.fits'
    fout = fb+'_nz.txt'
    common.mknz(fcd,fcr,fout,bs=dz_step,zmin=zmin,zmax=zmax)
    common.addnbar(fb,bs=dz_step,zmin=zmin,zmax=zmax,P0=P0,nran=nran,par=args.par,nproc=nproc,logger=logger)

def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    fn = Table(fitsio.read(flroot.replace('global','dvs_ro') +app))
    #c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    #gc = c.transform_to('galactic')
    sel_ngc = common.splitGC(fn)#gc.b > 0
    #outf_ngc = flroot.replace('/global/cfs/cdirs/desi/survey/catalogs/',os.getenv('SCRATCH')+'/')+'NGC_'+app
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS_scratchcp(fn[sel_ngc],outf_ngc,logger=logger)
    #outf_sgc = flroot.replace('/global/cfs/cdirs/desi/survey/catalogs/',os.getenv('SCRATCH')+'/')+'SGC_'+app
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS_scratchcp(fn[~sel_ngc],outf_sgc,logger=logger)

inds = np.arange(rm,rx)
if args.splitGC == 'y':
    #dirout = dirout.replace('/global/cfs/cdirs/desi/survey/catalogs/',os.getenv('SCRATCH')+'/')

    #if not os.path.exists(dirout):
    #    os.makedirs(dirout)
    #print('made '+dirout)
    splitGC(fb+'_','.dat')
    def _spran(rann):
        splitGC(fb+'_','.ran',rann)
    #inds = np.arange(nran)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            res = pool.map(_spran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _spran(rn)

if args.resamp == 'y':
    regions = ['NGC','SGC']        
    for reg in regions:
        flin = fb + '_'+reg    
        def _parfun_re(rannum):
            ct.clusran_resamp(flin,rannum,rcols=rcols)#,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)
        
        
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool() as pool:
                res = pool.map(_parfun_re, inds)
        else:
            for rn in range(rm,rx):
                _parfun_re(rn)
        fcr = flin+'_0_clustering.ran.fits'
        fcd = flin+'_clustering.dat.fits'
        fout = flin+'_nz.txt'
        common.mknz(fcd,fcr,fout,bs=dz_step,zmin=zmin,zmax=zmax)
        common.addnbar(flin,bs=dz_step,zmin=zmin,zmax=zmax,P0=P0,nran=nran,par=args.par,logger=logger)


if args.add_weight_ntile == 'y':
    bo = common.add_weight_ntile(fb, ranmin=rm, nran=rx, par=args.par)

   
