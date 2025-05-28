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
parser.add_argument("--mask_ran_nopriority", help="apply vetos to randoms except for priority; to be used with mocks",default='n')

parser.add_argument("--use_map_veto", help="string to include in full file name denoting whether map veto was applied",default='_HPmapcut')
parser.add_argument("--add_tlcomp", help="add completeness FRAC_TLOBS_TILES to randoms",default='n')



parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='')

parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--relax_zbounds", help="whether or not to use less restricted redshift bounds",default='y')
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
parser.add_argument("--replace_syscol",help="whether to replace any existing weight_sys with new",default='n')
parser.add_argument("--imsys_clus_fb_ran",help="add linear weight fits in fine redshift bins to randoms",default='n')



parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=10,type=int)

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

mainp = main(args.type,args.verspec,survey=args.survey,relax_zbounds=args.relax_zbounds)
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
if args.type == 'LGE':
    tarver = '3.0.0'
    progl = 'dark1b'
tardir = '/global/cfs/cdirs/desi/target/catalogs/dr9/'+tarver+'/targets/main/resolve/'
tarf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+type +'targetsDR9v'+tarver.strip('.')+'.fits'

mktar = True
if os.path.isfile(tarf) and redotar == False or len(type.split('-'))>1:
    mktar = False
#if type == 'BGS_BRIGHT':
#    mktar = False    

if mktar: #concatenate target files for given type, with column selection hardcoded
    import LSS.imaging.select_samples as ss
    ss.gather_targets(type,tardir,tarf,tarver,'main',progl,keys=keys)

mketar = True
etardir = '/global/cfs/cdirs/desi/survey/catalogs/extra_target_data/'+tarver+'/'
etarf = maindir+type +'targets_pixelDR9v'+tarver.strip('.')+'.fits'        
if os.path.isfile(etarf) and redotar == False: 
    mketar = False

if args.survey != 'main':
    mketar = False

if type == 'BGS_BRIGHT':
    mketar = False

if mketar: #concatenate target files for given type, with column selection hardcoded
    import LSS.imaging.select_samples as ss
    ss.gather_targets(type,etardir,etarf,tarver,'main',progl)

maxp = 3400
if type[:3] == 'LRG' or notqso == 'notqso':
    maxp = 3200
if type[:3] == 'BGS':
    maxp = 2100

       
if mkfulld:
    logf.write('creating full data catalogs for '+tp+' '+str(datetime.now()))
    logger.info('creating full data catalogs for '+tp)
    azf=''
    azfm = 'cumul'        
    if args.survey == 'DA02':
        #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+progl+'-cumulative.fits')
        #zmtlf = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/datcomb_'+progl+'_zmtl_zdone.fits')
        if type[:3] == 'ELG':
            azf = mainp.elgzf
            #azfm = 'hp'
        if type[:3] == 'QSO':
            azf = mainp.qsozf
            azfm = 'hp'
            if specrel == 'newQSOtemp_tagged':
                azf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/newQSOtemp_tagged/QSO_catalog.fits'
                azfm = 'cumul'
        dz = ldirspec+'datcomb_'+progl+'_tarspecwdup_zdone.fits' #new
        tlf = ldirspec+'Alltiles_'+progl+'_tilelocs.dat.fits'

    else:
        tracer_ts = type
        if type[:3] == 'ELG':
            tracer_ts = 'ELG'
        if type[:3] == 'BGS':
            tracer_ts = 'BGS_ANY'
        dz = ldirspec+'datcomb_'+tracer_ts+'_tarspecwdup_zdone.fits'
        tlf = None
        if type[:3] == 'ELG':
            azf = mainp.elgzf
        if type[:3] == 'QSO':
            azf = mainp.qsozf
    #if args.survey == 'main':        
    #    tlf = ldirspec+type+'_tilelocs.dat.fits'

 
    ftar = fitsio.read(tarf)   

    from desitarget import targetmask
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]
        desitarg='DESI_TARGET'
    
    maskcoll = False
    if args.survey != 'main':
        maskcoll = True
    ct.mkfulldat(dz,imbits,ftar,type,bit,dirout+type+notqso+'_full_noveto.dat.fits',tlf,survey=args.survey,maxp=maxp,azf=azf,azfm=azfm,desitarg=desitarg,specver=specrel,notqso=notqso,min_tsnr2=tsnrcut,badfib=mainp.badfib_td,badfib_status=mainp.badfib_status,mask_coll=maskcoll,logger=logger)


if args.add_veto == 'y':
    logf.write('added veto columns to data catalogs for '+tp+' '+str(datetime.now()))
    fin = dirout+type+notqso+'_full_noveto.dat.fits'
    common.add_veto_col(fin,ran=False,tracer_mask=type[:3].lower(),redo=True)#,rann=0
    for rn in range(rm,rx):
        fin = dirout+progl+'_'+str(rn)+'_full_noveto.ran.fits'
        common.add_veto_col(fin,ran=True,tracer_mask=type[:3].lower(),rann=rn)
        
if args.join_etar == 'y':
    logf.write('added extra target columns to data catalogs for '+tp+' '+str(datetime.now()))
    fin = dirout+type+notqso+'_full_noveto.dat.fits'
    common.join_etar(fin,type)

new_cols=mainp.new_cols#['STARDENS','HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z', 'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15', 'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK']
fid_cols=mainp.fid_cols#['EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
allmapcols = new_cols+fid_cols

if args.fillran == 'y':
    logf.write('filled randoms with imaging properties for '+tp+' '+str(datetime.now()))
    print('filling randoms with imaging properties')
    for ii in range(rm,rx):
        fn = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        #ct.addcol_ran(fn,ii)
        common.add_map_cols(fn,ii,new_cols=new_cols,fid_cols=fid_cols)
        print('done with '+str(ii))


if args.apply_veto == 'y':
    print('applying vetos')
    logf.write('applied vetos to data catalogs for '+tp+' '+str(datetime.now()))

    if args.ranonly != 'y':
        fin = dirout.replace('global','dvs_ro')+type+notqso+'_full_noveto.dat.fits'
        fout = dirout+type+notqso+'_full.dat.fits'
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp,reccircmasks=mainp.reccircmasks,logger=logger)
    print('data veto done, now doing randoms')
    def _parfun(rn):
        #fin = dirout.replace('global','dvs_ro')+type+notqso+'_'+str(rn)+'_full_noveto.ran.fits'
        fin = dirout.replace('global','dvs_ro')+progl+'_'+str(rn)+'_full_noveto.ran.fits'
        fout = dirout+type+notqso+'_'+str(rn)+'_full.ran.fits'
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp,reccircmasks=mainp.reccircmasks,logger=logger)
        print('random veto '+str(rn)+' done')
    if args.par == 'n':
        for rn in range(rm,rx):
            _parfun(rn)
    else:
        inds = np.arange(rm,rx)
        from multiprocessing import Pool
        (rx-rm)*2
        nproc = 9 #try this so doesn't run out of memory
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun, inds)


wzm = ''
if ccut is not None:
    wzm += ccut #you could change this to however you want the file names to turn out


tracer_clus = type+notqso+wzm

regl = ['_N','_S']    

lssmapdirout = dirout+'/hpmaps/'
nside = 256
if args.mkHPmaps == 'y':
    from LSS.imaging.sky_maps import create_pixweight_file, create_pixweight_file_allinone ,rancat_names_to_pixweight_name
    #logf.write('made healpix property maps for '+tp+' '+str(datetime.now()))
    if not os.path.exists(lssmapdirout):
        os.mkdir(lssmapdirout)
        common.printlog('made '+lssmapdirout,logger)
    lssmapdir = '/global/cfs/cdirs/desi/survey/catalogs/external_input_maps/'
    rancatname = dirout+tracer_clus+'_*_full.ran.fits'
    rancatlist = sorted(glob.glob(rancatname))
    #print(dirout)
    #print(rancatlist)
    fieldslist = allmapcols
    masklist = list(np.zeros(len(fieldslist),dtype=int))
    
    for reg in ['N','S']:
        outfn = lssmapdirout+tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
        create_pixweight_file(rancatlist, fieldslist, masklist, nside_out=nside,
                          lssmapdir=lssmapdir, outfn=outfn,reg=reg)    
    #regl = ['N','S']
    #outfn = lssmapdirout+tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_'
    #create_pixweight_file_allinone(rancatlist, fieldslist, masklist, nside_out=nside,
    #                      lssmapdir=lssmapdir, outfn=outfn,regl=regl)
if args.apply_map_veto == 'y':
    import healpy as hp
    tracer_clushp = tracer_clus
    #BGS_ANY and BGS_BRIGHT should essentially have same footprint
    if tracer_clus == 'BGS_ANY':
        tracer_clushp = 'BGS_BRIGHT'
    if 'ELG' in tracer_clus:
        tracer_clushp = 'ELG_LOPnotqso'
    mapn = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
    mapcuts = mainp.mapcuts
    common.printlog('apply map vetos',logger)
    if args.ranonly != 'y':
        fout = dirout+type+notqso+'_full.dat.fits'
        fin = fout.replace('global','dvs_ro')  
        fout = fout.replace('_full','_full_HPmapcut')      
        common.apply_map_veto(fin,fout,mapn,maps,mapcuts,logger=logger)
        common.printlog('data veto done, now doing randoms',logger)
    else:
        common.printlog('not doing data',logger)
    def _parfun(rn):
        fout = dirout+type+notqso+'_'+str(rn)+'_full.ran.fits'
        fin = fin = fout.replace('global','dvs_ro')   
        fout = fout.replace('_full','_full_HPmapcut')          
        common.apply_map_veto(fin,fout,mapn,maps,mapcuts,logger=logger)
        common.printlog('random veto '+str(rn)+' done',logger)
    if args.par == 'n':
        for rn in range(rm,rx):
            _parfun(rn)
    else:
        inds = np.arange(rm,rx)
        from multiprocessing import Pool
        
        nproc = 9 #try this so doesn't run out of memory
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun, inds)

if args.mask_ran_nopriority == 'y':
    import healpy as hp
    tracer_clushp = tracer_clus
    #BGS_ANY and BGS_BRIGHT should essentially have same footprint
    if tracer_clus == 'BGS_ANY':
        tracer_clushp = 'BGS_BRIGHT'
    if 'ELG' in tracer_clus:
        tracer_clushp = 'ELG_LOPnotqso'
    mapn = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
    mapcuts = mainp.mapcuts

    ran_fname_base = dirout+type+notqso+'_'
    def _mk_inputran(rann):
        outfn = ran_fname_base.replace('dvs_ro','global')+str(rann)+'_full_noPriveto_HPmapcut.ran.fits'
        if args.survey == 'Y1':
            infn = ran_fname_base+str(rann)+'_full_noveto.ran.fits'
        else:
            infn = dirout+prog.lower()+'_'+str(rann)+'_full_noveto.ran.fits'
        maxp = 10000 #we don't want to apply any priority cut
    
        masked_dat = common.apply_veto(infn,outfn,ebits=ebits,zmask=False,maxp=maxp,logger=logger,reccircmasks=mainp.reccircmasks,wo='n')
        masked_dat = common.apply_map_veto_arrays(masked_dat,mapn,maps,mapcuts)
        common.write_LSS_scratchcp(masked_dat,outfn,logger=logger)
    inds = np.arange(rm,rx)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool(processes=9) as pool:
            res = pool.map(_mk_inputran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _mk_inputran(rn)


if args.add_tlcomp == 'y':
    fl = dirout+tp+notqso+'_'
    def _parfun(ii):
        ct.add_tlobs_ran(fl,ii,args.use_map_veto)
    if args.par == 'n':
        for rn in range(rm,rx):
            _parfun(rn)
    else:
        inds = np.arange(rm,rx)
        from multiprocessing import Pool
        (rx-rm)*2
        nproc = 9 #try this so doesn't run out of memory
        with Pool(processes=nproc) as pool:
            res = pool.map(_parfun, inds)

rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL','WEIGHT_SN','WEIGHT_RF','TARGETID_DATA']#,'WEIGHT_FKP']#,'WEIGHT_RF']
if type[:3] == 'BGS':
    fcols = ['G','R','Z','W1','W2']
    for col in fcols:
        rcols.append('flux_'+col.lower()+'_dered')


if args.add_fs == 'y':
    fscols=['TARGETID','ABSMAG01_SDSS_G','ABSMAG01_SDSS_R']
    fsver = 'v2.0'
    fsrel = 'dr1'
    fsspecver = args.verspec
    logf.write('adding columns from fastspecfit version ' +fsver+' '+fsrel+' '+str(datetime.now()))
    if 'global' in dirout:
        #diro = copy(dirout)
        inroot = '/dvs_ro/cfs/cdirs/'
        outroot = '/global/cfs/cdirs/'
        infn = dirout.replace(outroot,'')+type+notqso+'_full'+args.use_map_veto+'.dat.fits'
        print(dirout)
    else:
        inroot = ''
        outroot = ''
        infn = dirout+type+notqso+'_full'+args.use_map_veto+'.dat.fits'
    common.join_with_fastspec(infn,fscols,inroot=inroot,\
    outroot=outroot,fsver=fsver,fsrel=fsrel,specrel=fsspecver,prog=progl)


if args.add_ke == 'y':
    logf.write('adding k+e columns '+tp+' '+str(datetime.now()))
    reglke = regl
    if args.survey != 'DA02':
        reglke = ['']
    kecols = ['REST_GMR_0P1','KCORR_R0P1','KCORR_G0P1','KCORR_R0P0','KCORR_G0P0','REST_GMR_0P0','EQ_ALL_0P0'\
    ,'EQ_ALL_0P1','REST_GMR_0P1','ABSMAG_RP0','ABSMAG_RP1'] 
    for col in kecols:
        rcols.append(col)

    for reg in reglke:
        fb = dirout+tracer_clus+reg
        if args.survey == 'DA02':
            fn = fb+'_clustering.dat.fits'
            zcol = 'Z'
        else:
            fn = fb+'_full.dat.fits'
            zcol = 'Z_not4clus'
        dat = Table(fitsio.read(fn))
        dat = common.add_dered_flux(dat,fcols)
        n_processes = 100
        from multiprocessing import Pool
        chunk_size = (len(dat)+n_processes)//n_processes
        list = []
        for i in range(0,n_processes):
            mini = i*chunk_size
            maxi = mini+chunk_size
            if maxi > len(dat):
                maxi = len(dat)
            list.append(dat[mini:maxi])
        
        def _wrapper(N):
            mini = N*chunk_size
            maxi = mini+chunk_size
            if maxi > len(dat):
                maxi = len(dat)
            idx = np.arange(mini,maxi)
            data = list[N]#Table()
            data['idx'] = idx
            #list[N] = common.add_ke(data,zcol='Z_not4clus')
            data = common.add_ke(data,zcol=zcol)
            return data

        with Pool(processes=n_processes+1) as pool:
            res = pool.map(_wrapper, np.arange(n_processes))
            #pool.map(_wrapper, np.arange(n_processes))

        res = vstack(res)#vstack(list)#
        res.sort('idx')
        res.remove_column('idx')
        print(len(res),len(dat))

        #if args.test == 'y':
        #    dat = dat[:10]
        #cols = list(dat.dtype.names)
        #if 'REST_GMR_0P1' in cols:
        #    print('appears columns are already in '+fn)
        #else:
        #    dat = common.add_ke(dat,zcol='Z_not4clus')
            #if args.test == 'n':
        common.write_LSS(res,fn,comments=['added k+e corrections'])

if 'BGS_BRIGHT-' in type:
#type == 'BGS_BRIGHT-21.5':# and args.survey == 'Y1': #and args.clusd == 'y':
    abmagcut = -float(type.split('-')[1])
    common.printlog('using ab mag cut '+str(abmagcut),logger)
    ffull = dirout+type+notqso+'_full'+args.use_map_veto+'.dat.fits'
    if os.path.isfile(ffull) == False or args.redoBGS215 == 'y':
        logf.write('making BGS_BRIGHT'+str(abmagcut)+' full data catalog for '+str(datetime.now()))
        fin = fitsio.read(dirout+'BGS_BRIGHT_full'+args.use_map_veto+'.dat.fits')
        if args.absmagmd == 'phot':
            sel = fin['ABSMAG_RP1'] < abmagcut
        if args.absmagmd == 'spec':
            sel = (fin['ABSMAG01_SDSS_R'] +0.97*fin['Z_not4clus']-.095) <abmagcut
            #sys.exit('need to code up using fastspecfit for abs mag selection!')
        if args.absmagmd == 'nok':
            #don't use any k-correction at all, yields ~constant density
            from LSS.tabulated_cosmo import TabulatedDESI
            cosmo = TabulatedDESI()
            dis_dc = cosmo.comoving_radial_distance
            z2use = np.copy(fin['Z_not4clus'])
            selz = z2use <= 0
            selz |= z2use > 2
            z2use[selz] = 2
            dm = 5.*np.log10(dis_dc(z2use)*(1.+z2use)) + 25.
            cfluxr = fin['FLUX_R']/fin['MW_TRANSMISSION_R']
            r_dered = 22.5 - 2.5*np.log10(cfluxr)
            abr = r_dered -dm
            if abmagcut == -21.5:
                sel = abr < -21.6 +0.15*z2use
            else:
                sel = abr < abmagcut
            sel &= z2use < 2
        common.write_LSS(fin[sel],ffull)

    
if args.add_bitweight == 'y':
    print('USE add_DR1_bitweights_all.py script instead!')
    #logf.write('added bitweights to data catalogs for '+tp+' '+str(datetime.now()))
    #fn = dirout+type+notqso+'_full'+args.use_map_veto+'.dat.fits'
    #print(fn)
    #ff = Table(fitsio.read(fn))
    #try:
    #    ff.remove_columns(['BITWEIGHTS_1','PROB_OBS_1','BITWEIGHTS_2','PROB_OBS_2'])
    #    print('removed ','BITWEIGHTS_1','PROBOBS_1','BITWEIGHTS_2','PROBOBS_2')
    #except:
    #    print('not removing 1/2 bitweights')
    #try:
    #    ff.remove_columns(['BITWEIGHTS','PROB_OBS'])
    #    print('removed ','BITWEIGHTS','PROB_OBS')
    #except:
    #    print('not removing bitweights')
    #
    #if type[:3] != 'BGS':
    #    bitf = fitsio.read(mainp.darkbitweightfile)
    #else:
    #    bitf = fitsio.read(mainp.brightbitweightfile)
    #ff = join(ff,bitf,keys=['TARGETID'],join_type='left')
    #common.write_LSS(ff,fn)#,comments='Added alt MTL info')

if args.swap20211212 == 'y':
    dirspec = '/global/cfs/cdirs/desi/spectro/redux/reproc_20211212_iron/tiles/cumulative/'
    tllist = [10376, 10380,  21386, 22541, 23406,  23414,  24518,  24567,2526,25266,26075,2840,2842,5642, 7207,  8621]
    
    datal = []
    keep_cols = ['TARGETID','Z','ZWARN','DELTACHI2']
    rename_cols = ['Z_not4clus','ZWARN','DELTACHI2']
    if type[:3] == 'ELG':
        eml = []
        keep_cols = ['TARGETID','Z','ZWARN','DELTACHI2','OII_FLUX','OII_FLUX_IVAR']
        rename_cols = ['Z_not4clus','ZWARN','DELTACHI2','OII_FLUX','OII_FLUX_IVAR']

    for tl in tllist:
        rr = fitsio.read(dirspec+str(tl)+'/20211212/redrock-9-'+str(tl)+'-thru20211212.fits')
        datal.append(rr)
        if type[:3] == 'ELG':
            em = fitsio.read(dirspec+str(tl)+'/20211212/emline-9-'+str(tl)+'-thru20211212.fits')
            eml.append(em)
    data = np.concatenate(datal)
    data = Table(data)
    if type[:3] == 'ELG':
        emd = np.concatenate(eml)
        c2add = ['OII_FLUX','OII_FLUX_IVAR']
        for col in c2add:
            data[col] = emd[col]
    data.keep_columns(keep_cols)
    data['Z'].name = 'Z_not4clus'
    for col in rename_cols:
        data[col].name = col+'_4swap'
    swap_cols = ['Z_not4clus_4swap','ZWARN_4swap','DELTACHI2_4swap']
    if type[:3] == 'ELG':
        o2c = np.log10(data['OII_FLUX_4swap'] * np.sqrt(data['OII_FLUX_IVAR_4swap']))+0.2*np.log10(data['DELTACHI2_4swap'])
        data['o2c_4swap'] = o2c
        swap_cols = ['Z_not4clus_4swap','ZWARN_4swap','DELTACHI2_4swap','OII_FLUX_4swap','OII_FLUX_IVAR_4swap','o2c_4swap']
    fn = dirout+type+notqso+'_full'+args.use_map_veto+'.dat.fits'
    
    ff = fitsio.read(fn)
    ff = join(ff,data,keys=['TARGETID'],join_type='left')
    for col in rename_cols:
        ff[col+'_4swap'] = ff[col+'_4swap'].filled(999999)

    sel = ff['Z_not4clus_4swap'] != 999999
    sel &= np.isin(ff['TILEID'],tllist)
    sel &= ff['ZWARN'] != 999999 #don't flip any original not observed
    print('rows to have redshift values replaced:')
    print(len(ff[sel]))
    for col in swap_cols:
        ff[col.replace('_4swap','')][sel] = ff[col][sel]
    ff.remove_columns(swap_cols)
    common.write_LSS(ff,fn)

if args.fixzwarn == 'y':
    fn = dirout+type+notqso+'_full'+args.use_map_veto+'.dat.fits'    
    ff = Table(fitsio.read(fn))
    sel = ff['FIBER'] == 999999
    sel &= ff['ZWARN'] != 999999
    ff['ZWARN'][sel] = 999999
    sel = ff['FIBER'] == 999999
    sela = ff['ZWARN'] == 999999
    print(len(ff[sel]),len(ff[sela]))
    common.write_LSS(ff,fn)
    
#if mkclusran and mkclusdat:
#    print('doing clustering randoms')
#    for ii in range(rm,rx):
#        ct.mkclusran(dirin+type+notqso+'_',dirout+tracer_clus+'_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits)#,ntilecut=ntile,ccut=ccut)


if args.add_weight_zfail == 'y':
    logf.write('regressing and adding WEIGHT_ZFAIL to data catalogs for '+tp+' '+str(datetime.now()))
    readpars = False
    if args.readpars == 'y':
        readpars = True
    if type[:3] == 'QSO':
        ct.add_zfail_weight2fullQSO(ldirspec,version,mainp.qsozf,tsnrcut=tsnrcut,readpars=readpars,logger=logger) 
    else:
        ct.add_zfail_weight2full(dirout,tp=type+notqso,tsnrcut=tsnrcut,readpars=readpars,logger=logger)   

if args.usemaps == None:
    fit_maps = mainp.fit_maps
else:
    fit_maps = [mapn for mapn in args.usemaps]

zl = (zmin,zmax)
#fit_maps = ['EBV_CHIANG_SFDcorr','STARDENS','HALPHA','EBV_MPF_Mean_FW15','BETA_ML','HI','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z']



tpstr = tracer_clus
if 'BGS_BRIGHT' in tracer_clus:
    tpstr = 'BGS_BRIGHT'
nside = 256


if type[:3] == 'ELG':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.1),(1.1,1.6)]
    else:
        zrl = [(0.8,1.6)]
if type[:3] == 'QSO':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.3),(1.3,2.1),(2.1,3.5)] 
    else:
        zrl = [(0.8,3.5)]   
if type[:3] == 'LRG':
    if args.imsys_zbin == 'y':
        #zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
        zrl = [(0.3,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2)] 
    else:
        zrl = [(0.4,1.1)]
    zsysmin = 0.4
    zsysmax = 1.1
    if args.relax_zbounds == 'y':
        zsysmax = 1.2      
if 'BGS_BRIGHT-' in type:
    zrl = [(0.1,0.4)]
elif type[:3] == 'BGS':
    zrl = [(0.01,0.5)]
    #zmin = 0.01
    #zmax = 0.5    



if args.prepsysnet == 'y' or args.regressis == 'y' or args.imsys == 'y' or args.imsys_clus == 'y':
    
    debv = common.get_debv()
    
    sky_g,sky_r,sky_z = common.get_skyres()
    #ebvn_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/test/initial_corrected_ebv_map_nside_64.fits'
    #ebvn = fitsio.read(ebvn_fn)
    #debv = ebvn['EBV_NEW'] - ebvn['EBV_SFD']
    #debv64 = make_hp(debv, ebvn['HPXPIXEL'], nside=64, fill_with=hp.UNSEEN)
    #debv256 = hp.ud_grade(debv64, 256)
    #debv256_nest = hp.reorder(debv256,r2n=True)
    #debv = Table()
    #debv['EBV_DIFFRZ'] = debv256_nest

if args.imsys == 'y':
    from LSS.imaging import densvar
    #regl = ['_DN','_DS','','_N','_S']
    #wzm = ''
    #fit_maps = ['STARDENS','EBV','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
    
       
    #rcols.append('WEIGHT_SYSEB')   
    fname = os.path.join(dirout, tracer_clus+'_full'+args.use_map_veto+'.dat.fits')
    dat = Table(fitsio.read(fname))
    selgood = common.goodz_infull(tp[:3],dat)
    selobs = dat['ZWARN'] != 999999
    dat = dat[selgood&selobs]
    ranl = []
    for i in range(0,args.nran4imsys):#int(args.maxr)):
        ran = fitsio.read(os.path.join(dirout, tpstr+'_'+str(i)+'_full'+args.use_map_veto+'.ran.fits'), columns=['RA', 'DEC','PHOTSYS']) 
        ranl.append(ran)
    rands = np.concatenate(ranl)
    common.printlog('combined randoms',logger)
    syscol = 'WEIGHT_IMLIN'
    regl = ['S','N']
    if args.type == 'QSO':
        regl = ['DES','SnotDES','N']

    dat[syscol] = np.ones(len(dat))
    for reg in regl:
        regu = reg
        if 'DES' in reg:
            regu = 'S'
        pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'_'+regu+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            #if 'EBV_DIFF_'+ec in fit_maps: 
            sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        if 'EBV_DIFF_MPF' in fit_maps:
            sys_tab['EBV_DIFF_MPF'] = sys_tab['EBV'] - sys_tab['EBV_MPF_Mean_FW15']
        #for bnd in ['G','R','Z']:
        if 'SKY_RES_G' in fit_maps:
            sys_tab['SKY_RES_G'] = sky_g[regu]
        if 'SKY_RES_R' in fit_maps:
            sys_tab['SKY_RES_R'] = sky_r[regu]
        if 'SKY_RES_Z' in fit_maps:
            sys_tab['SKY_RES_Z'] = sky_z[regu]
           
        #seld = dat['PHOTSYS'] == reg
        #selr = rands['PHOTSYS'] == reg

        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            #fb = dirout+tracer_clus+reg
            #fcr = fb+'_0_clustering.ran.fits'
            #rd = fitsio.read(fcr)
            #fcd = fb+'_clustering.dat.fits'
            #dd = Table.read(fcd)
            
            common.printlog('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax),logger)
            if type == 'LRG' and args.usemaps is None:
                fitmapsbin = mainp.fit_maps_allebv
                #if reg == 'N':
                #    fitmapsbin = fit_maps
                #else:
                #    if zmax == 0.6:
                #        fitmapsbin = mainp.fit_maps46s
                #    if zmax == 0.8:
                #        fitmapsbin = mainp.fit_maps68s
                #    if zmax == 1.1:
                #        fitmapsbin = mainp.fit_maps81s
            else:
                fitmapsbin = fit_maps
            use_maps = fitmapsbin
            wsysl = densvar.get_imweight(dat,rands,zmin,zmax,reg,fitmapsbin,use_maps,sys_tab=sys_tab,zcol='Z_not4clus',figname=dirout+tracer_clus+'_'+reg+'_'+str(zmin)+str(zmax)+'_linimsysfit.png')
            sel = wsysl != 1
            dat[syscol][sel] = wsysl[sel]
            #dd['WEIGHT'][sel] *= wsysl[sel]
            #dd.write(fcd,overwrite=True,format='fits')
    common.write_LSS_scratchcp(dat,fname,logger=logger)

if args.prepsysnet == 'y':
    logf.write('preparing data to run sysnet regression for '+tp+' '+str(datetime.now())+'\n')
    if not os.path.exists(dirout+'/sysnet'):
        os.mkdir(dirout+'/sysnet')
        print('made '+dirout+'/sysnet')    

    from LSS.imaging import sysnet_tools
    
    #allsky_rands = None
    #if args.use_allsky_rands == 'y':
    #    print('using randoms allsky for frac_area')
        #ranl = []
        #randir = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
        #for i in range(0,18):
        #    logger.info("reading allsky randoms "+str(i))
        #    ran = fitsio.read(randir+f'randoms-allsky-1-{i}.fits',columns=['RA','DEC','PHOTSYS'])
        #    ranl.append(ran)
        #allsky_rands = np.concatenate(ranl)
        
    
    #_HPmapcut'
    dat = fitsio.read(os.path.join(dirout, tracer_clus+'_full'+args.use_map_veto+'.dat.fits'))
    import time
    time_start = time.time()
    print('Loading randoms...')
    ranl = []
    if args.par == 'y':
        def read_catalogs(randoms_path):
            return fitsio.read(randoms_path,columns=['RA', 'DEC','PHOTSYS'])
        rands_paths=[]
        for i in range(0,18):
            rands_paths.append(os.path.join(dirout.replace('global','dvs_ro'), tpstr+'_'+str(i)+'_full'+args.use_map_veto+'.ran.fits')) 
        from multiprocessing import Pool
        with Pool() as pool:
            ranl = pool.map(read_catalogs, rands_paths)
    else:
        for i in range(0,18):
            ran = fitsio.read(os.path.join(dirout.replace('global','dvs_ro'), tpstr+'_'+str(i)+'_full'+args.use_map_veto+'.ran.fits'), columns=['RA', 'DEC','PHOTSYS']) 
            ranl.append(ran)
    rands = np.concatenate(ranl)
    print('Randoms loaded and stacked', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
    regl = ['N','S']
    
    for zl in zrl:
        zw = ''
        zmin,zmax=zl[0],zl[1]
        if args.imsys_zbin == 'y':
            zw = str(zmin)+'_'+str(zmax)
        for reg in regl:
            if type == 'LRG':
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
            tpmap = tpstr
            if 'ELG' in tpstr:
                tpmap = 'ELG_LOPnotqso'
            pwf = lssmapdirout+tpmap+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
            sys_tab = Table.read(pwf)
            cols = list(sys_tab.dtype.names)
            for col in cols:
                if 'DEPTH' in col:
                    bnd = col.split('_')[-1]
                    sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
            for ec in ['GR','RZ']:
                if 'EBV_DIFF_'+ec in fit_maps: 
                    sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
            if 'EBV_DIFF_MPF' in fit_maps:
                sys_tab['EBV_DIFF_MPF'] = sys_tab['EBV'] - sys_tab['EBV_MPF_Mean_FW15']

            seld = dat['PHOTSYS'] == reg
            selr = rands['PHOTSYS'] == reg
            if args.use_allsky_rands == 'y':
                allsky_fn = f"/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/allsky_rpix_{reg}_nran18_nside256_ring.fits"
                allsky_rands = fitsio.read(allsky_fn)
                allrands = allsky_rands['RANDS_HPIX'] # randoms count per hp pixel
            #    selr_all = allsky_rands['PHOTSYS'] == reg
            #    allrands = allsky_rands[selr_all]
            else:
                allrands = None
            common.printlog(f"{tpstr} {reg} z{zmin}-{zmax}: {fitmapsbin}",logger)
            prep_table = sysnet_tools.prep4sysnet(dat[seld], rands[selr], sys_tab, zcolumn='Z_not4clus', allsky_rands=allrands, 
                                                  zmin=zl[0], zmax=zl[1], nran_exp=None, nside=nside, nest=True, use_obiwan=False,
                                                  columns=fitmapsbin,wtmd='fracz',tp=args.type[:3])
            fnout = dirout+'/sysnet/prep_'+tracer_clus+zw+'_'+reg+'.fits'
            common.write_LSS(prep_table,fnout)

if args.regressis == 'y':
    logf.write('adding regressis weights to data catalogs for '+tp+' '+str(datetime.now())+'\n')
    #from regressis, must be installed
    from regressis import DESIFootprint,DR9Footprint

    from LSS.imaging import regressis_tools as rt
    dirreg = dirout+'/regressis_data'
    
#     if type[:3] == 'ELG':
#         zl = (0.8,1.5)
#     if type[:3] == 'QSO':
#         zl = (0.8,2.1)#,(2.1,3.5)]    
#     if type[:3] == 'LRG':
#         zl = (0.4,1.1)
#     if type[:3] == 'BGS':
#         zl = (0.1,0.5)  


    if not os.path.exists(dirreg):
        os.mkdir(dirreg)
        print('made '+dirreg)   
    #pwf = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits'   
    sgf = '/global/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_'+str(nside)+'.npy' 
    dr9_footprint = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    desi_footprint = DESIFootprint(nside)
    suffix_tracer = ''
    suffix_regressor = ''

    param = dict()
    param['data_dir'] = dirreg
    param['output_dir'] = dirreg
    param['use_median'] = False
    param['use_new_norm'] = False
    if tracer_clus[:3] == 'QSO':
        param['regions'] = ['North', 'South', 'Des']
    else:
        param['regions'] = ['North', 'South_mid_ngc', 'South_mid_sgc']
    max_plot_cart = 300

    cut_fracarea = False
    seed = 42
    #fit_maps = None
    '''
    Map choices are:
    'EBV_CHIANG_SFDcorr','STARDENS','HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z',
    'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15',
    'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1',
    'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK'
    'EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1',
    'PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z'
    '''
    #if tracer_clus[:3] == 'BGS':# or tracer_clus[:3] == 'ELG':
    #    fit_maps = ['STARDENS','EBV','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
        #fit_maps = ['STARDENS','HI','BETA_ML','PSFDEPTH_G', 'PSFDEPTH_R','PSFDEPTH_Z','PSFDEPTH_W1','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']

    logf.write('using fit maps '+str(fit_maps)+'\n')
    feature_names_ext=None
    pw_out_fn_root = dirout+'/regressis_data/'+tracer_clus+'feature_data_'
    regl = ['N','S']
    
    for reg in regl:
        pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        pw_out_fn = pw_out_fn_root+reg+'.fits'
    
        print(pw_out_fn)
        sys_tab.write(pw_out_fn,overwrite=True,format='fits')

    #pixweight_data = Table.read(pwf)
    #if 'EBV_DIFFRZ' in fit_maps: 
    #    pixweight_data['EBV_DIFFRZ'] = debv['EBV_DIFFRZ']
    #for ec in ['GR','RZ']:
    #    if 'EBV_DIFF_'+ec in fit_maps: 
    #        pixweight_data['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        
    use_sgr=False
    if 'SGR' in fit_maps:
        use_sgr = True
        fit_maps.remove('SGR')

    for zl in zrl:    
        zw = str(zl[0])+'_'+str(zl[1])
        #rt.save_desi_data_full(dirout, 'main', tracer_clus, nside, dirreg, zl,foot=dr9_footprint,nran=18)
    
        common.printlog('computing regressis weight in mode '+args.regmode+'for '+tracer_clus+zw,logger)
        #logf.write('computing RF regressis weight for '+tracer_clus+zw+'\n')
        
        rt.get_desi_data_full_compute_weight(dirout, 'main', tracer_clus, nside, dirreg, zl, param,foot=dr9_footprint,nran=18,\
        suffix_tracer=suffix_tracer, suffix_regressor=suffix_regressor, cut_fracarea=cut_fracarea, seed=seed,\
         max_plot_cart=max_plot_cart,pixweight_path=pw_out_fn_root,pixmap_external=debv,sgr_stream_path=sgf,\
         feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext,use_map_veto=args.use_map_veto,regressor=args.regmode)
        #rt._compute_weight('main', tracer_clus+zw, dr9_footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, max_plot_cart,pixweight_path=pw_out_fn,pixmap_external=debv,sgr_stream_path=sgf,feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext)

if args.add_regressis == 'y':
    from LSS.imaging import densvar
    from regressis import PhotoWeight
    fb = dirout+tracer_clus
    fcd = fb+'_full'+args.use_map_veto+'.dat.fits'
    dd = Table.read(fcd)
    dd['WEIGHT_'+args.regmode] = np.ones(len(dd))

    for zl in zrl:    
        print(zl)
        zw = str(zl[0])+'_'+str(zl[1])

        fnreg = dirout+'/regressis_data/main_'+tracer_clus+zw+'_256/'+args.regmode+'/main_'+tracer_clus+zw+'_imaging_weight_256.npy'
        #fracarea = np.load(dirout+'/regressis_data/main_'+tracer_clus+zw+'_fracarea_256.npy')
        #selarea = fracarea*0 == 0
        #rfw = np.load(fnreg,allow_pickle=True)
        #rfpw = rfw.item()['map']
        #maskreg = rfw.item()['mask_region']
        #regl_reg = list(maskreg.keys())
        #for reg in regl_reg:
        #    mr = maskreg[reg]
        #    norm = np.mean(rfpw[mr&selarea])
        #    print(reg,norm)
        #    rfpw[mr] /= norm
        rfpw = PhotoWeight.load(fnreg)
        #print(np.mean(rfpw))
        #dth,dphi = densvar.radec2thphi(dd['RA'],dd['DEC'])
        #dpix = densvar.hp.ang2pix(densvar.nside,dth,dphi,nest=densvar.nest)
        #drfw = rfpw[dpix]
        
        selz = dd['Z_not4clus'] > zl[0]
        selz &= dd['Z_not4clus'] <= zl[1]
        dd['WEIGHT_'+args.regmode][selz] = rfpw(dd['RA'][selz], dd['DEC'][selz], normalize_map=True)#drfw[selz]
        #norm = 
        print(np.mean(dd['WEIGHT_'+args.regmode][selz]))
    #comments = []
    #comments.append("Using regressis for WEIGHT_SYS")
    logf.write('added RF regressis weight for '+tracer_clus+zw+'\n')

    common.write_LSS(dd,fcd)#,comments)

if args.add_regressis_ext == 'y':
    if tracer_clus != 'QSO':
        sys.exit('only QSO supported for using weights Edmond calculated!')
    from LSS.imaging import densvar
    from regressis import PhotoWeight
    fb = dirout+tracer_clus
    fcd = fb+'_full.dat.fits'
    dd = Table.read(fcd)
    dd['WEIGHT_SYS'] = np.ones(len(dd))
    
    wsys_low_z = PhotoWeight.load('/global/homes/e/edmondc/CFS/Imaging_weight/Y1/Y1_QSO_low_z_128/RF/Y1_QSO_low_z_imaging_weight_128.npy')
    wsys_high_z = PhotoWeight.load('/global/homes/e/edmondc/CFS/Imaging_weight/Y1/Y1_QSO_high_z_128/RF/Y1_QSO_high_z_imaging_weight_128.npy')

    sel_low_z = dd['Z_not4clus'] <= 1.3
    dd['WEIGHT_SYS'][sel_low_z] = wsys_low_z(dd['RA'][sel_low_z], dd['DEC'][sel_low_z], normalize_map=True)
    print(np.mean(dd['WEIGHT_SYS'][sel_low_z]))
    dd['WEIGHT_SYS'][~sel_low_z] = wsys_high_z(dd['RA'][~sel_low_z], dd['DEC'][~sel_low_z], normalize_map=True)
    print(np.mean(dd['WEIGHT_SYS'][~sel_low_z]))
    logf.write('added RF regressis weights that Edmond calculated for '+tracer_clus+'\n')    
    common.write_LSS(dd,fcd)
 

if args.add_sysnet == 'y':
    logf.write('adding sysnet weights to data catalogs for '+tp+' '+str(datetime.now())+'\n')
    from LSS.imaging import densvar
    import healpy as hp
    fn_full = dirout+tracer_clus+'_full'+args.use_map_veto+'.dat.fits'
    dd = Table.read(fn_full)
    dd['WEIGHT_SN'] = np.ones(len(dd))
    dth,dphi = densvar.radec2thphi(dd['RA'],dd['DEC'])
    dpix = hp.ang2pix(256,dth,dphi)

    regl_sysnet = ['N','S']
    for reg in regl_sysnet:
        for zl in zrl:
            zw = ''
            if args.imsys_zbin == 'y':
                zw = str(zl[0])+'_'+str(zl[1])
            sn_weights = fitsio.read(dirout+'/sysnet/'+tracer_clus+zw+'_'+reg+'/nn-weights.fits')
            pred_counts = np.mean(sn_weights['weight'],axis=1)
            #pix_weight = np.mean(pred_counts)/pred_counts
            #pix_weight = np.clip(pix_weight,0.5,2.)
            pix_weight = 1./pred_counts
            pix_weight = pix_weight / pix_weight.mean()
            pix_weight = np.clip(pix_weight,0.5,2.)
            sn_pix = sn_weights['hpix']
            hpmap = np.ones(12*256*256)
            for pix,wt in zip(sn_pix,pix_weight):
                hpmap[pix] = wt
        
            sel = dd['PHOTSYS'] == reg
            selz = dd['Z_not4clus'] > zl[0]
            selz &= dd['Z_not4clus'] <= zl[1]

            #print(np.sum(sel))
            dd['WEIGHT_SN'][sel&selz] = hpmap[dpix[sel&selz]]
            if tracer_clus == 'ELG_LOPnotqso':
                if zl[0] == 0.8:
                    selz = dd['Z_not4clus'] <= zl[0]
                if zl[1] == 1.6:
                    selz = dd['Z_not4clus'] > zl[1]
                dd['WEIGHT_SN'][sel&selz] = hpmap[dpix[sel&selz]]
        #assign weights to galaxies outside the z ranges
        if tracer_clus == 'ELG_LOPnotqso':
            zwl = '0.8_1.1'
            
    #print(np.min(dd['WEIGHT_SYS']),np.max(dd['WEIGHT_SYS']),np.std(dd['WEIGHT_SYS']))
    comments = []
    comments.append("Using sysnet for WEIGHT_SYS")

    common.write_LSS(dd,fn_full,comments)
        
    


utlid = False
if args.ran_utlid == 'y':
    utlid = True


#needs to happen before randoms so randoms can get z and weights
weightileloc=True
if args.compmd == 'altmtl':
    weightileloc = False
if mkclusdat:
    ct.mkclusdat(dirout+type+notqso,weightileloc,tp=type,dchi2=dchi2,zmin=zmin,zmax=zmax,correct_zcmb='n',wsyscol=args.imsys_colname,use_map_veto=args.use_map_veto,extradir=args.extra_clus_dir)#,ntilecut=ntile,ccut=ccut)

nzcompmd = 'ran'
if args.compmd == 'altmtl':
    nzcompmd = args.compmd


inds = np.arange(rm,rx)
if mkclusran:
    print('doing clustering randoms (possibly a 2nd time to get sys columns in)')
#     tsnrcol = 'TSNR2_ELG'
#     tsnrcut = 0
#    
#     if type[:3] == 'ELG':
#         #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
#         tsnrcut = 80
#     if type == 'LRG':
#         #dchi2 = 16  
#         tsnrcut = 80  
#     if type[:3] == 'BGS':
#         tsnrcol = 'TSNR2_BGS'
#         dchi2 = 40
#         tsnrcut = 1000
    ranin = dirin + args.type + notqso + '_'
    if 'BGS_BRIGHT' in args.type:
        ranin = dirin + 'BGS_BRIGHT' + notqso + '_'

    clus_arrays = [fitsio.read(dirout +args.extra_clus_dir+ type + notqso+'_clustering.dat.fits')]
    def _parfun_cr(ii):
        ct.mkclusran(ranin,dirout+tracer_clus+'_',ii,rcols=rcols,ebits=ebits,utlid=utlid,clus_arrays=clus_arrays,use_map_veto=args.use_map_veto,compmd=nzcompmd,logger=logger,extradir=args.extra_clus_dir,tp=type)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_parfun_cr, inds)

    else:
        for ii in inds:#range(rm,rx):
            _parfun_cr(ii)
        #,ntilecut=ntile,ccut=ccut)

if args.NStoGC == 'y':
    fb = dirout+tracer_clus+'_'
    ct.clusNStoGC(fb, rx - rm)#,par=args.par)


if type == 'QSO':
    #zmin = 0.6
    #zmax = 4.5
    dz = 0.02
    P0 = 6000
    
else:    
    dz = 0.01
    #zmin = 0.01
    #zmax = 1.61

if type[:3] == 'LRG':
    P0 = 10000
if type[:3] == 'ELG':
    P0 = 4000
if type[:3] == 'BGS':
    P0 = 7000

nran = rx-rm
regions = ['NGC', 'SGC']

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


if args.splitGC == 'y':
    fb = dirout+args.extra_clus_dir+tracer_clus+'_'
   # ct.splitclusGC(fb, args.maxr - args.minr,par=args.par)   
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


if args.resamp == 'y':
            
    for reg in regions:
        flin = dirout + tracer_clus + '_'+reg    
        def _parfun(rannum):
            ct.clusran_resamp(flin,rannum,rcols=rcols)#,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)
        
        
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool() as pool:
                res = pool.map(_parfun, inds)
        else:
            for rn in range(rm,rx):
                _parfun(rn)
    
#allreg = ['N','S','NGC', 'SGC']
#allreg = ['NGC','SGC']
if args.nz == 'y':
    for reg in regions:#allreg:
        fb = dirout+args.extra_clus_dir+tracer_clus+'_'+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.txt'
        zmu = (zmin//dz)*dz #make sure some integer multiple of dz
        common.printlog('minimum z used in n(z) calc is '+str(zmu),logger)
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmu,zmax=zmax,compmd=nzcompmd)
        common.addnbar(fb,bs=dz,zmin=zmu,zmax=zmax,P0=P0,nran=nran,par=args.par,compmd=nzcompmd)

if args.addNtileweight2full == 'y':
    froot = dirout+tracer_clus
    if args.survey == 'Y1':
        nproc = 9
    if args.survey == 'DA2':
        nproc = 9
    common.add_weight_ntile(froot,logger=logger,ranmin=rm,nran=rx,par=args.par,extradir=args.extra_clus_dir,tp=type,nproc=nproc)

if args.imsys_clus == 'y':
    from LSS.imaging import densvar
    
       
    #rcols.append('WEIGHT_SYSEB')   
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname))
    dat = vstack([dat_sgc,dat_ngc])
    foutname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_clustering.dat.fits')
    ranl = []
    for i in range(0,args.nran4imsys):#int(args.maxr)):
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)

    rands = np.concatenate(ranl)
    syscol = 'WEIGHT_IMLIN_CLUS'
    regl = ['S','N']
    if args.type == 'QSO':
        regl = ['DES','SnotDES','N']
    dat[syscol] = np.ones(len(dat))
    for reg in regl:
        regu = reg
        if reg == 'DES' or reg == 'SnotDES':
            regu = 'S'
        pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'_'+regu+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        #seld = dat['PHOTSYS'] == reg
        selr = rands['PHOTSYS'] == reg

        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            
            print('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax))
            if type == 'LRG':
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
            wsysl = densvar.get_imweight(dat,rands,zmin,zmax,reg,fitmapsbin,use_maps,sys_tab=sys_tab,zcol='Z',modoutname = dirout+tracer_clus+'_'+reg+'_'+str(zmin)+str(zmax)+'_linfitparam.txt',figname=dirout+tracer_clus+'_'+reg+'_'+str(zmin)+str(zmax)+'_linclusimsysfit.png',wtmd='clus')
            sel = wsysl != 1
            dat[syscol][sel] = wsysl[sel]
    #attach data to NGC/SGC catalogs, write those out
    dat.keep_columns(['TARGETID',syscol])
    if syscol in dat_ngc.colnames:
        dat_ngc.remove_column(syscol)
    dat_ngc = join(dat_ngc,dat,keys=['TARGETID'])
    common.write_LSS_scratchcp(dat_ngc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits'),logger=logger)
    if syscol in dat_sgc.colnames:
        dat_sgc.remove_column(syscol)
    dat_sgc = join(dat_sgc,dat,keys=['TARGETID'])
    common.write_LSS_scratchcp(dat_sgc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits'),logger=logger)

if args.imsys_clus_ran == 'y':
    #do randoms
    syscol = 'WEIGHT_IMLIN_CLUS'
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    dat = vstack([dat_sgc,dat_ngc])
    dat.rename_column('TARGETID','TARGETID_DATA')
    regl = ['NGC','SGC']
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_'+reg+'_'+str(rann)+'_clustering.ran.fits')
            ran = Table(fitsio.read(ran_fn))
            if syscol in ran.colnames:
                ran.remove_column(syscol)
            ran = join(ran,dat,keys=['TARGETID_DATA'])
            common.write_LSS_scratchcp(ran,ran_fn,logger=logger)

    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _add2ran(rn)
            
if args.imsys_clus_fb == 'y':
    from LSS.imaging import densvar
    
       
    #rcols.append('WEIGHT_SYSEB')   
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname))
    dat = vstack([dat_sgc,dat_ngc])
    foutname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_clustering.dat.fits')
    ranl = []
    for i in range(0,args.nran4imsys):#int(args.maxr)):
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)

    rands = np.concatenate(ranl)
    syscol = 'WEIGHT_IMLIN_CLUS'
    regl = ['S','N']
    if args.type == 'QSO':
        regl = ['DES','SnotDES','N']
    dat[syscol] = np.ones(len(dat))
    for reg in regl:
        regu = reg
        if reg == 'DES' or reg == 'SnotDES':
            regu = 'S'
        pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'_'+regu+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        #seld = dat['PHOTSYS'] == reg
        selr = rands['PHOTSYS'] == reg
        dz = 0.1
        zm = zsysmin
        zx = zm + dz
        fitmapsbin = fit_maps
        while zm < zsysmax:
            zx = zm + dz
            zx = round(zx,1)
            print('getting weights for region '+reg+' and '+str(zm)+'<z<'+str(zx))
            if type == 'LRG':
                fitmapsbin = mainp.fit_maps_all
            else:
                fitmapsbin = fit_maps
            use_maps = fitmapsbin
            wsysl = densvar.get_imweight(dat,rands,zm,zx,reg,fitmapsbin,use_maps,sys_tab=sys_tab,zcol='Z',modoutname = dirout+args.extra_clus_dir+tracer_clus+'_'+reg+'_'+str(zm)+str(zx)+'_linfitparam.txt',figname=dirout+tracer_clus+'_'+reg+'_'+str(zm)+str(zx)+'_linclusimsysfit.png',wtmd='clus')
            sel = wsysl != 1
            dat[syscol][sel] = wsysl[sel]

            zm = zx
            #zm += dz
            #zm = round(zm,1)
    #attach data to NGC/SGC catalogs, write those out
    dat.keep_columns(['TARGETID',syscol])
    if syscol in dat_ngc.colnames:
        dat_ngc.remove_column(syscol)
    dat_ngc = join(dat_ngc,dat,keys=['TARGETID'])
    if args.replace_syscol == 'y':
        dat_ngc['WEIGHT'] /= dat_ngc['WEIGHT_SYS']
        dat_ngc['WEIGHT_SYS'] = dat_ngc[syscol]
        dat_ngc['WEIGHT'] *= dat_ngc['WEIGHT_SYS']
    common.write_LSS_scratchcp(dat_ngc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits'),logger=logger)
    if syscol in dat_sgc.colnames:
        dat_sgc.remove_column(syscol)
    dat_sgc = join(dat_sgc,dat,keys=['TARGETID'])
    if args.replace_syscol == 'y':
        dat_sgc['WEIGHT'] /= dat_sgc['WEIGHT_SYS']
        dat_sgc['WEIGHT_SYS'] = dat_sgc[syscol]
        dat_sgc['WEIGHT'] *= dat_sgc['WEIGHT_SYS']

    common.write_LSS_scratchcp(dat_sgc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits'),logger=logger)

if args.imsys_clus_fb_ran == 'y':
    #do randoms
    syscol = 'WEIGHT_IMLIN_CLUS'
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    dat = vstack([dat_sgc,dat_ngc])
    dat.rename_column('TARGETID','TARGETID_DATA')
    regl = ['NGC','SGC']
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_'+reg+'_'+str(rann)+'_clustering.ran.fits')
            ran = Table(fitsio.read(ran_fn))
            if syscol in ran.colnames:
                ran.remove_column(syscol)
            ran = join(ran,dat,keys=['TARGETID_DATA'])
            if args.replace_syscol == 'y':
                ran['WEIGHT'] /= ran['WEIGHT_SYS']
                ran['WEIGHT_SYS'] = ran[syscol]
                ran['WEIGHT'] *= ran['WEIGHT_SYS']
            
            common.write_LSS_scratchcp(ran,ran_fn,logger=logger)

    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _add2ran(rn)


#if args.nz == 'y':
    
#     if zmask:
#         wzm = 'zmask_'
#     if rcut is not None:
#         wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
#     if ntile > 0:
#         wzm += '_ntileg'+str(ntilecut)+'_'    

#    regl = ['_DN','_DS','','_N','_S']
    
    
#    for reg in regl:
#        fb = dirout+tracer_clus+reg
#        fcr = fb+'_0_clustering.ran.fits'
#        fcd = fb+'_clustering.dat.fits'
#        fout = fb+'_nz.txt'
#        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
#        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0)

# if args.FKPfull == 'y':
#     
#     fb = dirout+tracer_clus
#     fbr = fb
#     if type == 'BGS_BRIGHT-21.5':
#         fbr = dirout+'BGS_BRIGHT'
# 
#     fcr = fbr+'_0_full.ran.fits'
#     fcd = fb+'_full.dat.fits'
#     nz = common.mknz_full(fcd,fcr,type[:3],bs=dz,zmin=zmin,zmax=zmax)
#     common.addFKPfull(fcd,nz,type[:3],bs=dz,zmin=zmin,zmax=zmax,P0=P0)
# 
# if args.nzfull == 'y':
#     fb = dirout+tracer_clus
#     fbr = fb
#     if type == 'BGS_BRIGHT-21.5':
#         fbr = dirout+'BGS_BRIGHT'
#     fcr = fbr+'_0_full.ran.fits'
#     fcd = fb+'_full.dat.fits'
#     zmax = 1.6
#     zmin = 0.01
#     bs = 0.01
#     if type[:3] == 'QSO':
#         zmax = 4
#         bs = 0.02
#     for reg in regl:
#         reg = reg.strip('_')
#         common.mknz_full(fcd,fcr,type[:3],bs,zmin,zmax,randens=2500.,write='y',reg=reg)    
#         nzf = np.loadtxt(fb+'_full_'+reg+'_nz.txt').transpose()
#         plt.plot(nzf[0],nzf[3],label=reg)
#     plt.xlabel('redshift')
#     plt.ylabel('n(z) (h/Mpc)^3')
#     plt.legend()
#     plt.grid()
#     if tracer_clus == 'ELG_LOPnotqso':
#         plt.ylim(0,0.001)
#     if tracer_clus == 'BGS_BRIGHT':
#         plt.yscale('log')
#         plt.xlim(0,0.6)
#         plt.ylim(1e-5,0.15)
#     if tracer_clus == 'BGS_BRIGHT-21.5':
#         plt.xlim(0,0.5)
#     plt.title(tracer_clus)
#     plt.savefig(dirout+'plots/'+tracer_clus+'_nz.png')

# if args.addnbar_ran == 'y':
#     utlid_sw = ''
#     if utlid:
#         utlid_sw = '_utlid'
#     
#     for reg in regl:
#         fb = dirout+tracer_clus+utlid_sw+reg
#         common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,add_data=False,ran_sw=utlid_sw)
# 
# 
# if args.swapz == 'y':
#     import LSS.blinding_tools as blind
#     for reg in regl:
#         fb = dirout+tracer_clus+reg+'_clustering.dat.fits'
#         data = Table(fitsio.read(fb))
#         blind.swap_z(data,fb,frac=0.01)        
