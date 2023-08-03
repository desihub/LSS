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
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--basedir_blind", help="base directory for output for blinded catalogs, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--redotar", help="remake the target file for the particular type (needed if, e.g., the requested columns are changed)",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--add_veto", help="add veto column for given type, matching to targets",default='n')
parser.add_argument("--join_etar", help="whether or not to join to the target files with extra brick pixel info",default='n')
parser.add_argument("--apply_veto", help="apply vetos for imaging, priorities, and hardware failures",default='n')
parser.add_argument("--mkHPmaps", help="make healpix maps for imaging properties using sample randoms",default='n')
parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--nzfull", help="get n(z) from full files",default='n')

parser.add_argument("--FKPfull", help="add FKP weights to full catalogs",default='n')
parser.add_argument("--addnbar_ran", help="just add nbar/fkp to randoms",default='n')
parser.add_argument("--add_ke", help="add k+e corrections for BGS data to clustering catalogs",default='n')
parser.add_argument("--add_fs", help="add rest frame info from fastspecfit",default='n')
parser.add_argument("--absmagmd", help="whether to use purely photometry+z based rest frame info or fastspecfit",choices=['spec','phot'],default='spec')

parser.add_argument("--blinded", help="are we running on the blinded full catalogs?",default='n')

parser.add_argument("--prepsysnet",help="prepare data to get sysnet weights for imaging systematics?",default='n')
parser.add_argument("--add_sysnet",help="add sysnet weights for imaging systematics to full files?",default='n')
parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='n')



parser.add_argument("--regressis",help="RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis",help="add RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis_ext",help="add RF weights for imaging systematics, calculated elsewhere",default='n')
parser.add_argument("--imsys_nside",help="healpix nside used for imaging systematic regressions",default=256,type=int)


parser.add_argument("--add_weight_zfail",help="add weights for redshift systematics to full file?",default='n')
parser.add_argument("--add_bitweight",help="add info from the alt mtl",default='n')


parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

parser.add_argument("--par", help="run different random number in parallel?",default='y')

parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)
parser.add_argument("--ranonly",help="if y, only operate on randoms when applying vetos",default='n')
parser.add_argument("--readpars",help="set to y for certain steps if you want to read from previous fits",default='n')


#options not typically wanted
parser.add_argument("--imsys",help="add weights for imaging systematics using eboss method?",default='n')
parser.add_argument("--ran_utlid", help="cut randoms so that they only have 1 entry per tilelocid",default='n')
parser.add_argument("--swapz", help="if blinded, swap some fraction of redshifts?",default='n')


args = parser.parse_args()
print(args)

type = args.type
tp = type
basedir = args.basedir
version = args.version
specrel = args.verspec
ntile = args.ntile
ccut = args.ccut
rm = int(args.minr)
rx = int(args.maxr)


print('running catalogs for tracer type '+type)

redotar = False
if args.redotar == 'y':
    redotar = True

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
        
    
if mkfulld:
    print('making "full" catalog file for data')    
    
    
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
    print('making clustering catalog for data')
    
if args.clusran == 'y':
    mkclusran = True
    
if mkclusran:
    print('making clustering catalog for randoms, files '+str(rm)+ ' through '+str(rx))
    print('(if running all, consider doing in parallel)')  

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
    print('made '+ldirspec)
    
if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    print('made '+ldirspec+'LSScats')

dirout = ldirspec+'LSScats/'+version+'/'
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
    print('made '+dirout)    

tarver = '1.1.1'
tardir = '/global/cfs/cdirs/desi/target/catalogs/dr9/'+tarver+'/targets/main/resolve/'
tarf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+type +'targetsDR9v'+tarver.strip('.')+'.fits'

mktar = True
if os.path.isfile(tarf) and redotar == False or type == 'BGS_BRIGHT-21.5':
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
        dz = ldirspec+'datcomb_'+type+'_tarspecwdup_zdone.fits'
        tlf = ldirspec+type+'_tilelocs.dat.fits'
        if type[:3] == 'ELG':
            azf = mainp.elgzf
        if type[:3] == 'QSO':
            azf = mainp.qsozf
    if args.survey == 'Y1':
        tlf = None

 
    ftar = fitsio.read(tarf)   

    from desitarget import targetmask
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]
        desitarg='DESI_TARGET'
    
    maskcoll = False
    if args.survey == 'Y1':
        maskcoll = True
    ct.mkfulldat(dz,imbits,ftar,type,bit,dirout+type+notqso+'_full_noveto.dat.fits',tlf,survey=args.survey,maxp=maxp,azf=azf,azfm=azfm,desitarg=desitarg,specver=specrel,notqso=notqso,min_tsnr2=tsnrcut,badfib=mainp.badfib,mask_coll=maskcoll)

if args.add_bitweight == 'y':
    logf.write('added bitweights to data catalogs for '+tp+' '+str(datetime.now()))
    fn = dirout+type+notqso+'_full_noveto.dat.fits'
    ff = fitsio.read(fn)
    if type[:3] != 'BGS':
        bitf = fitsio.read(mainp.darkbitweightfile)
    else:
        bitf = fitsio.read(mainp.brightbitweightfile)
    ff = join(ff,bitf,keys=['TARGETID'],join_type='left')
    common.write_LSS(ff,fn,comments='Added alt MTL info')


if args.add_veto == 'y':
    logf.write('added veto columns to data catalogs for '+tp+' '+str(datetime.now()))
    fin = dirout+type+notqso+'_full_noveto.dat.fits'
    common.add_veto_col(fin,ran=False,tracer_mask=type[:3].lower(),redo=True)#,rann=0
    for rn in range(rm,rx):
        fin = dirout+type+notqso+'_'+str(rn)+'_full_noveto.ran.fits'
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
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp,reccircmasks=mainp.reccircmasks)
    print('data veto done, now doing randoms')
    def _parfun(rn):
        fin = dirout.replace('global','dvs_ro')+type+notqso+'_'+str(rn)+'_full_noveto.ran.fits'
        fout = dirout+type+notqso+'_'+str(rn)+'_full.ran.fits'
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp,reccircmasks=mainp.reccircmasks)
        print('random veto '+str(rn)+' done')
    if args.par == 'n':
        for rn in range(rm,rx):
            _parfun(rn)
    else:
        inds = np.arange(rm,rx)
        from multiprocessing import Pool
        with Pool(processes=nran*2) as pool:
            res = pool.map(_parfun, inds)


wzm = ''
if ccut is not None:
    wzm += ccut #you could change this to however you want the file names to turn out

if type == 'BGS_BRIGHT-21.5' and args.survey == 'Y1':
    ffull = dirout+type+notqso+'_full.dat.fits'
    if os.path.isfile(ffull) == False:
        logf.write('making BGS_BRIGHT-21.5 full data catalog for '+str(datetime.now()))
        fin = fitsio.read(dirout+'BGS_BRIGHT_full.dat.fits')
        if args.absmagmd == 'phot':
            sel = fin['ABSMAG_RP1'] < -21.5
        if args.absmagmd == 'spec':
            sel = (fin['ABSMAG_SDSS_R'] +0.97*fin['Z_not4clus']-.095) < -21.5
            #sys.exit('need to code up using fastspecfit for abs mag selection!')
        common.write_LSS(fin[sel],ffull)

tracer_clus = type+notqso+wzm

regl = ['_N','_S']    
#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    ct.mkclusdat(dirout+type+notqso,tp=type,dchi2=dchi2,tsnrcut=tsnrcut,zmin=zmin,zmax=zmax,ccut=ccut)#,ntilecut=ntile)

lssmapdirout = dirout+'/hpmaps/'
if args.mkHPmaps == 'y':
    from LSS.imaging.sky_maps import create_pixweight_file, rancat_names_to_pixweight_name
    logf.write('made healpix property maps for '+tp+' '+str(datetime.now()))
    if not os.path.exists(lssmapdirout):
        os.mkdir(lssmapdirout)
        print('made '+lssmapdirout)
    lssmapdir = '/global/cfs/cdirs/desi/survey/catalogs/external_input_maps/'
    rancatname = dirout+tracer_clus+'_*_full.ran.fits'
    rancatlist = sorted(glob.glob(rancatname))
    fieldslist = allmapcols
    masklist = list(np.zeros(len(fieldslist),dtype=int))
    nside = 256
    for reg in ['N','S']:
        outfn = lssmapdirout+tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
        create_pixweight_file(rancatlist, fieldslist, masklist, nside_out=nside,
                          lssmapdir=lssmapdir, outfn=outfn,reg=reg)    

rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']#,'WEIGHT_FKP']#,'WEIGHT_RF']
if type[:3] == 'BGS':
    fcols = ['G','R','Z','W1','W2']
    for col in fcols:
        rcols.append('flux_'+col.lower()+'_dered')


if args.add_fs == 'y':
    fscols=['TARGETID','ABSMAG_SDSS_G','ABSMAG_SDSS_R']
    fsver = 'v1.0'
    fsrel = 'dr1'
    fsspecver = args.verspec
    logf.write('adding columns from fastspecfit version ' +fsver+' '+fsrel+' '+str(datetime.now()))
    if 'global' in dirout:
        #diro = copy(dirout)
        inroot = '/dvs_ro/cfs/cdirs/'
        outroot = '/global/cfs/cdirs/'
        infn = dirout.replace(outroot,'')+type+notqso+'_full.dat.fits'
        print(dirout)
    else:
        inroot = ''
        outroot = ''
        infn = dirout+type+notqso+'_full.dat.fits'
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
    



if mkclusran and mkclusdat:
    print('doing clustering randoms')
    for ii in range(rm,rx):
        ct.mkclusran(dirin+type+notqso+'_',dirout+tracer_clus+'_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits)#,ntilecut=ntile,ccut=ccut)


if args.add_weight_zfail == 'y':
    logf.write('regressing and adding WEIGHT_ZFAIL to data catalogs for '+tp+' '+str(datetime.now()))
    readpars = False
    if args.readpars == 'y':
        readpars = True
    if type[:3] == 'QSO':
        ct.add_zfail_weight2fullQSO(ldirspec,version,mainp.qsozf,tsnrcut=tsnrcut,readpars=readpars) 
    else:
        ct.add_zfail_weight2full(dirout,tp=type+notqso,tsnrcut=tsnrcut,readpars=readpars)   


if args.imsys == 'y':
    from LSS.imaging import densvar
    #regl = ['_DN','_DS','','_N','_S']
    #wzm = ''
    fit_maps = ['STARDENS','EBV','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
    use_maps = fit_maps
    if type[:3] == 'ELG':
        zrl = [(0.6,0.8),(0.8,1.1),(1.1,1.5)]
    if type[:3] == 'QSO':
        zrl = [(0.8,1.6),(1.6,2.1),(2.1,3.5)]    
    if type[:3] == 'LRG':
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]    
    if type[:3] == 'BGS':
        zrl = [(0.1,0.5)]    
       
    rcols.append('WEIGHT_SYSEB')   
    
    for reg in regl:
        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            fb = dirout+tracer_clus+reg
            fcr = fb+'_0_clustering.ran.fits'
            rd = fitsio.read(fcr)
            fcd = fb+'_clustering.dat.fits'
            dd = Table.read(fcd)
            dd['WEIGHT_SYSEB'] = np.ones(len(dd))
            print('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax))
            wsysl = densvar.get_imweight(dd,rd,zmin,zmax,fit_maps,use_maps,plotr=False)
            sel = wsysl != 1
            dd['WEIGHT_SYSEB'][sel] = wsysl[sel]
            #dd['WEIGHT'][sel] *= wsysl[sel]
            dd.write(fcd,overwrite=True,format='fits')

zl = (zmin,zmax)
#fit_maps = ['EBV_CHIANG_SFDcorr','STARDENS','HALPHA','EBV_MPF_Mean_FW15','BETA_ML','HI','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z']
fit_maps = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','EBV_DIFF_GR','EBV_DIFF_RZ','HI']
if tracer_clus[:3] == 'LRG':
    fit_maps.append('PSFDEPTH_W1')
#    fit_maps = ['STARDENS','HI','BETA_ML','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
if tracer_clus[:3] == 'QSO':
    fit_maps = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','EBV_DIFF_GR','EBV_DIFF_RZ']
    fit_maps.append('PSFDEPTH_W1')
    fit_maps.append('PSFDEPTH_W2')
    #fit_maps = ['EBV', 'STARDENS',
    #             'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
    #             'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2']

tpstr = tracer_clus
if tracer_clus == 'BGS_BRIGHT-21.5':
    tpstr = 'BGS_BRIGHT'
nside = 256
pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'.fits'

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
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
    else:
        zrl = [(0.4,1.1)]  
if type[:3] == 'BGS':
    zrl = [(0.1,0.4)]    


if args.prepsysnet == 'y' or args.regressis == 'y':
    
    def make_hp(value, hpix, nside, fill_with=np.nan):
        """ A Function to create a HEALPix map
        """
        m_ = np.zeros(12*nside*nside)
        m_[:] = fill_with
        m_[hpix] = value
    
        return m_

    import healpy as hp

    dirmap = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
    nside = 256#64
    nest = False
    eclrs = ['gr','rz']
    debv = Table()
    for ec in eclrs:
        ebvn = fitsio.read(dirmap+'v0_desi_ebv_'+ec+'_'+str(nside)+'.fits')
        debv_a = ebvn['EBV_DESI_'+ec.upper()]-ebvn['EBV_SFD']
        debv_a = hp.reorder(debv_a,r2n=True)
        debv['EBV_DIFF_'+ec.upper()] = debv_a

    #ebvn_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/test/initial_corrected_ebv_map_nside_64.fits'
    #ebvn = fitsio.read(ebvn_fn)
    #debv = ebvn['EBV_NEW'] - ebvn['EBV_SFD']
    #debv64 = make_hp(debv, ebvn['HPXPIXEL'], nside=64, fill_with=hp.UNSEEN)
    #debv256 = hp.ud_grade(debv64, 256)
    #debv256_nest = hp.reorder(debv256,r2n=True)
    #debv = Table()
    #debv['EBV_DIFFRZ'] = debv256_nest



if args.prepsysnet == 'y':
    logf.write('preparing data to run sysnet regression for '+tp+' '+str(datetime.now())+'\n')
    if not os.path.exists(dirout+'/sysnet'):
        os.mkdir(dirout+'/sysnet')
        print('made '+dirout+'/sysnet')    

    from LSS.imaging import sysnet_tools
    dat = fitsio.read(os.path.join(dirout, f'{tpstr}'+'_full.dat.fits'))
    ranl = []
    for i in range(0,18):
        ran = fitsio.read(os.path.join(dirout, f'{tpstr}'+'_'+str(i)+'_full.ran.fits'), columns=['RA', 'DEC','PHOTSYS']) 
        ranl.append(ran)
    rands = np.concatenate(ranl)
    regl = ['N','S']
    sys_tab = Table.read(pwf)
    #if 'EBV_DIFFRZ' in fit_maps:
    #    sys_tab['EBV_DIFFRZ'] = debv['EBV_DIFFRZ']
    for ec in ['GR','RZ']:
        if 'EBV_DIFF_'+ec in fit_maps: 
            sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]

    for zl in zrl:
        zw = ''
        if args.imsys_zbin == 'y':
            zw = str(zl[0])+'_'+str(zl[1])
        for reg in regl:
            seld = dat['PHOTSYS'] == reg
            selr = rands['PHOTSYS'] == reg
        
            prep_table = sysnet_tools.prep4sysnet(dat[seld], rands[selr], sys_tab, zcolumn='Z_not4clus', zmin=zl[0], zmax=zl[1], nran_exp=None,
                    nside=nside, nest=True, use_obiwan=False, columns=fit_maps,wtmd='fracz',tp=args.type[:3])
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
    pixweight_data = Table.read(pwf)
    #if 'EBV_DIFFRZ' in fit_maps: 
    #    pixweight_data['EBV_DIFFRZ'] = debv['EBV_DIFFRZ']
    for ec in ['GR','RZ']:
        if 'EBV_DIFF_'+ec in fit_maps: 
            pixweight_data['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        
    use_sgr=False
    if 'SGR' in fit_maps:
        use_sgr = True
        fit_maps.remove('SGR')
    pw_out_fn = dirout+'/regressis_data/'+tracer_clus+'feature_data.fits'
    print(pw_out_fn)
    pixweight_data.write(pw_out_fn,overwrite=True,format='fits')

    for zl in zrl:    
        zw = str(zl[0])+'_'+str(zl[1])
        #rt.save_desi_data_full(dirout, 'main', tracer_clus, nside, dirreg, zl,foot=dr9_footprint,nran=18)
    
        print('computing RF regressis weight for '+tracer_clus+zw)
        logf.write('computing RF regressis weight for '+tracer_clus+zw+'\n')
        rt.get_desi_data_full_compute_weight(dirout, 'main', tracer_clus, nside, dirreg, zl, param,foot=dr9_footprint,nran=18,\
        suffix_tracer=suffix_tracer, suffix_regressor=suffix_regressor, cut_fracarea=cut_fracarea, seed=seed,\
         max_plot_cart=max_plot_cart,pixweight_path=pw_out_fn,pixmap_external=debv,sgr_stream_path=sgf,\
         feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext)
        #rt._compute_weight('main', tracer_clus+zw, dr9_footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, max_plot_cart,pixweight_path=pw_out_fn,pixmap_external=debv,sgr_stream_path=sgf,feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext)

if args.add_regressis == 'y':
    from LSS.imaging import densvar
    from regressis import PhotoWeight
    fb = dirout+tracer_clus
    fcd = fb+'_full.dat.fits'
    dd = Table.read(fcd)
    dd['WEIGHT_SYS'] = np.ones(len(dd))

    for zl in zrl:    
        print(zl)
        zw = str(zl[0])+'_'+str(zl[1])

        fnreg = dirout+'/regressis_data/main_'+tracer_clus+zw+'_256/RF/main_'+tracer_clus+zw+'_imaging_weight_256.npy'
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
        dd['WEIGHT_SYS'][selz] = rfpw(dd['RA'][selz], dd['DEC'][selz], normalize_map=True)#drfw[selz]
        #norm = 
        print(np.mean(dd['WEIGHT_SYS'][selz]))
    comments = []
    comments.append("Using regressis for WEIGHT_SYS")
    logf.write('added RF regressis weight for '+tracer_clus+zw+'\n')

    common.write_LSS(dd,fcd,comments)

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
    fn_full = dirout+tracer_clus+'_full.dat.fits'
    dd = Table.read(fn_full)
    dd['WEIGHT_SYS'] = np.ones(len(dd))
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
            pix_weight = np.mean(pred_counts)/pred_counts
            pix_weight = np.clip(pix_weight,0.5,2.)
            sn_pix = sn_weights['hpix']
            hpmap = np.ones(12*256*256)
            for pix,wt in zip(sn_pix,pix_weight):
                hpmap[pix] = wt
        
            sel = dd['PHOTSYS'] == reg
            selz = dd['Z_not4clus'] > zl[0]
            selz &= dd['Z_not4clus'] <= zl[1]

            #print(np.sum(sel))
            dd['WEIGHT_SYS'][sel&selz] = hpmap[dpix[sel&selz]]
    #print(np.min(dd['WEIGHT_SYS']),np.max(dd['WEIGHT_SYS']),np.std(dd['WEIGHT_SYS']))
    comments = []
    comments.append("Using sysnet for WEIGHT_SYS")

    common.write_LSS(dd,fn_full,comments)
        
    


utlid = False
if args.ran_utlid == 'y':
    utlid = True
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

    for ii in range(rm,rx):
        ct.mkclusran(dirin+type+notqso+'_',dirout+tracer_clus+'_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits,utlid=utlid)#,ntilecut=ntile,ccut=ccut)

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
    

if args.nz == 'y':
    
#     if zmask:
#         wzm = 'zmask_'
#     if rcut is not None:
#         wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
#     if ntile > 0:
#         wzm += '_ntileg'+str(ntilecut)+'_'    

#    regl = ['_DN','_DS','','_N','_S']
    
    
    for reg in regl:
        fb = dirout+tracer_clus+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.txt'
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0)

if args.FKPfull == 'y':
    
    fb = dirout+tracer_clus
    fbr = fb
    if type == 'BGS_BRIGHT-21.5':
        fbr = dirout+'BGS_BRIGHT'

    fcr = fbr+'_0_full.ran.fits'
    fcd = fb+'_full.dat.fits'
    nz = common.mknz_full(fcd,fcr,type[:3],bs=dz,zmin=zmin,zmax=zmax)
    common.addFKPfull(fcd,nz,type[:3],bs=dz,zmin=zmin,zmax=zmax,P0=P0)

if args.nzfull == 'y':
    fb = dirout+tracer_clus
    fbr = fb
    if type == 'BGS_BRIGHT-21.5':
        fbr = dirout+'BGS_BRIGHT'
    fcr = fbr+'_0_full.ran.fits'
    fcd = fb+'_full.dat.fits'
    zmax = 1.6
    zmin = 0.01
    bs = 0.01
    if type[:3] == 'QSO':
        zmax = 4
        bs = 0.02
    for reg in regl:
        reg = reg.strip('_')
        common.mknz_full(fcd,fcr,type[:3],bs,zmin,zmax,randens=2500.,write='y',reg=reg)    
        nzf = np.loadtxt(fb+'_full_'+reg+'_nz.txt').transpose()
        plt.plot(nzf[0],nzf[3],label=reg)
    plt.xlabel('redshift')
    plt.ylabel('n(z) (h/Mpc)^3')
    plt.legend()
    plt.grid()
    if tracer_clus == 'ELG_LOPnotqso':
        plt.ylim(0,0.001)
    if tracer_clus == 'BGS_BRIGHT':
        plt.yscale('log')
        plt.xlim(0,0.6)
        plt.ylim(1e-5,0.15)
    if tracer_clus == 'BGS_BRIGHT-21.5':
        plt.xlim(0,0.5)
    plt.title(tracer_clus)
    plt.savefig(dirout+'plots/'+tracer_clus+'_nz.png')

if args.addnbar_ran == 'y':
    utlid_sw = ''
    if utlid:
        utlid_sw = '_utlid'
    
    for reg in regl:
        fb = dirout+tracer_clus+utlid_sw+reg
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,add_data=False,ran_sw=utlid_sw)


if args.swapz == 'y':
    import LSS.blinding_tools as blind
    for reg in regl:
        fb = dirout+tracer_clus+reg+'_clustering.dat.fits'
        data = Table(fitsio.read(fb))
        blind.swap_z(data,fb,frac=0.01)        