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
from desitarget import targetmask
#from regressis, must be installed
from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.imaging.select_samples as ss
from LSS.globals import main
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')
parser.add_argument("--redotar", help="remake the target file for the particular type (needed if, e.g., the requested columns are changed)",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--add_veto", help="add veto column for given type, matching to targets",default='n')
parser.add_argument("--apply_veto", help="apply vetos for imaging, priorities, and hardware failures",default='n')
parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--imsys",help="add weights for imaging systematics?",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--blinded", help="are we running on the blinded full catalogs?",default='n')
parser.add_argument("--swapz", help="if blinded, swap some fraction of redshifts?",default='n')

parser.add_argument("--regressis",help="RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis",help="add RF weights for imaging systematics?",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)


args = parser.parse_args()
print(args)

type = args.type
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

mainp = main(args.type,args.verspec)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

#wt = mtld['FAPRGRM'] == progl
#wt &= mtld['SURVEY'] == 'main'
#wt &= mtld['ZDONE'] == 'true'
#mtld = mtld[wt]
#print('there are '+str(len(mtld))+' tiles')


#columns to select from target sample
keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
       'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
       'MASKBITS','WISEMASK_W1','WISEMASK_W2', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET','BGS_TARGET','SHAPE_R']



#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)
    
if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    print('made '+ldirspec+'LSScats')

dirout = ldirspec+'LSScats/'+version+'/'
dirin = dirout
if args.blinded == 'y':
    
    dirout += 'blinded/'

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    

tarver = '1.1.1'
tardir = '/global/cfs/cdirs/desi/target/catalogs/dr9/'+tarver+'/targets/main/resolve/'
tarf = maindir+type +'targetsDR9v'+tarver.strip('.')+'.fits'

mktar = True
if os.path.isfile(tarf) and redotar == False:
    mktar = False
#if type == 'BGS_BRIGHT':
#    mktar = False    

if mktar: #concatenate target files for given type, with column selection hardcoded
    ss.gather_targets(type,tardir,maindir,tarver,'main',progl,keys=keys)
        
        
if mkfulld:
    azf=''
    azfm = 'cumul'        
    if specrel == 'daily':
        dz = ldirspec+'datcomb_'+type+'_tarspecwdup_zdone.fits'
        tlf = ldirspec+type+'_tilelocs.dat.fits'
        if type[:3] == 'ELG':
            #azf = '/global/cfs/cdirs/desi/users/raichoor/spectro/daily/main-elg-daily-tiles-cumulative.fits'
            azf = ldirspec+'emlin_catalog.fits'
        if type[:3] == 'QSO':
            azf =ldirspec+'QSO_catalog.fits'
    #if specrel == 'daily':
        #specf = Table.read(ldirspec+'datcomb_'+progl+'_spec_zdone.fits')
    else:
        #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+progl+'-cumulative.fits')
        #zmtlf = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/datcomb_'+progl+'_zmtl_zdone.fits')
        if type[:3] == 'ELG':
            azf = mainp.elgzf
            azfm = 'hp'
        if type[:3] == 'QSO':
            azf = mainp.qsozf
            azfm = 'hp'
        dz = ldirspec+'datcomb_'+progl+'_tarspecwdup_zdone.fits' #new
        tlf = ldirspec+'Alltiles_'+progl+'_tilelocs.dat.fits'


 
    ftar = fitsio.read(tarf)   

    
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]
        desitarg='DESI_TARGET'
    
    ct.mkfulldat(dz,imbits,ftar,type,bit,dirout+type+notqso+'zdone_full_noveto.dat.fits',tlf,azf=azf,azfm=azfm,desitarg=desitarg,specver=specrel,notqso=notqso)

if args.add_veto == 'y':
    fin = dirout+type+notqso+'zdone_full_noveto.dat.fits'
    common.add_veto_col(fin,ran=False,tracer_mask=type[:3].lower(),redo=True)#,rann=0
    for rn in range(rm,rx):
        fin = dirout+type+notqso+'zdone_'+str(rn)+'_full_noveto.ran.fits'
        common.add_veto_col(fin,ran=True,tracer_mask=type[:3].lower(),rann=rn)
        
if args.apply_veto == 'y':
    print('applying vetos')
    maxp = 3400
    if type[:3] == 'LRG' or notqso == 'notqso':
        maxp = 3200
    if type[:3] == 'BGS':
        maxp = 2100
    fin = dirout+type+notqso+'zdone_full_noveto.dat.fits'
    fout = dirout+type+notqso+'zdone_full.dat.fits'
    common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp)
    print('data veto done, now doing randoms')
    for rn in range(rm,rx):
        fin = dirout+type+notqso+'zdone_'+str(rn)+'_full_noveto.ran.fits'
        fout = dirout+type+notqso+'zdone_'+str(rn)+'_full.ran.fits'
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp)
        print('random veto '+str(rn)+' done')

dchi2 = 9
tsnrcut = 0
if type[:3] == 'ELG':
	dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
	tsnrcut = 80
	zmin = 0.8
	zmax = 1.6
if type == 'LRG':
	dchi2 = 16  
	tsnrcut = 80
	zmin = 0.4
	zmax = 1.1  
if type[:3] == 'BGS':
	dchi2 = 40
	tsnrcut = 1000
	zmin = 0.1
	zmax = 0.5
if type == 'QSO':
	zmin = 0.8
	zmax = 3.5
        

regl = ['_N','_S']    
#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    ct.mkclusdat(dirout+type+notqso+'zdone_',tp=type,dchi2=dchi2,tsnrcut=tsnrcut,zmin=zmin,zmax=zmax)#,ntilecut=ntile,ccut=ccut)

if args.fillran == 'y':
    print('filling randoms with imaging properties')
    for ii in range(rm,rx):
        fn = dirout+type+notqso+'zdone_'+str(ii)+'_full.ran.fits'
        ct.addcol_ran(fn,ii)
        print('done with '+str(ii))

if mkclusran:
    print('doing clustering randoms')
    tsnrcol = 'TSNR2_ELG'
    tsnrcut = 0
   
    if type[:3] == 'ELG':
        #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
        tsnrcut = 80
    if type == 'LRG':
        #dchi2 = 16  
        tsnrcut = 80  
    if type[:3] == 'BGS':
        tsnrcol = 'TSNR2_BGS'
        dchi2 = 40
        tsnrcut = 1000
    rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']#,'WEIGHT_FKP']#,'WEIGHT_RF'
    if type[:3] == 'BGS':
        rcols.append('flux_r_dered')

    for ii in range(rm,rx):
        ct.mkclusran(dirin+type+notqso+'zdone_',dirout+type+notqso+'zdone_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits)#,ntilecut=ntile,ccut=ccut)


if args.imsys == 'y':
    from LSS.imaging import densvar
    #regl = ['_DN','_DS','','_N','_S']
    wzm = ''
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
       
        
    
    for reg in regl:
        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            fb = dirout+type+notqso+'zdone'+wzm+reg
            fcr = fb+'_0_clustering.ran.fits'
            rd = fitsio.read(fcr)
            fcd = fb+'_clustering.dat.fits'
            dd = Table.read(fcd)
            print('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax))
            wsysl = densvar.get_imweight(dd,rd,zmin,zmax,fit_maps,use_maps,plotr=False)
            sel = wsysl != 1
            dd['WEIGHT_SYS'][sel] = wsysl[sel]
            dd['WEIGHT'][sel] *= wsysl[sel]
            dd.write(fcd,overwrite=True,format='fits')

zl = (zmin,zmax)
if args.regressis == 'y':
    from LSS.imaging import regressis_tools as rt
    dirreg = dirout+'/regressis_data'
    nside = 256
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
    pwf = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits'    
    rt.save_desi_data(dirout, 'main', type+notqso, nside, dirreg, zl,regl=regl) 
    dr9_footprint = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=True, cut_desi=False)

    suffix_tracer = ''
    suffix_regressor = ''

    param = dict()
    param['data_dir'] = dirreg
    param['output_dir'] = dirreg
    param['use_median'] = False
    param['use_new_norm'] = False
    param['regions'] = ['North', 'South', 'Des']
    max_plot_cart = 1000

    cut_fracarea = False
    seed = 42

    rt._compute_weight('main', type+notqso, dr9_footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, param, max_plot_cart,pixweight_path=pwf)

if args.add_regressis == 'y':
    from LSS.imaging import densvar
    fnreg = dirout+'/regressis_data/main_'+type+notqso+'_256/RF/main_'+type+notqso+'_imaging_weight_256.npy'
    rfw = np.load(fnreg,allow_pickle=True)
    rfpw = rfw.item()['map']
    #regl = ['_DN','_DS','','_N','_S']
    for reg in regl:
        fb = dirout+type+notqso+'zdone'+reg
        fcd = fb+'_clustering.dat.fits'
        dd = Table.read(fcd)
        dth,dphi = densvar.radec2thphi(dd['RA'],dd['DEC'])
        dpix = densvar.hp.ang2pix(densvar.nside,dth,dphi,nest=densvar.nest)
        drfw = rfpw[dpix]
        dd['WEIGHT_SYS'] = drfw
        dd['WEIGHT'] *= dd['WEIGHT_SYS']
        #dd.write(fcd,format='fits',overwrite=True)
        comments = ["DA02 'clustering' LSS catalog for data, "+reg+" entries are only for data with good redshifts with "+str(zmin)+'<z<'+str(zmax)]
        comments = ["Using regressis for WEIGHT_SYS"]

        common.write_LSS(dd,fcd,comments)

    
    
rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']#,'WEIGHT_FKP']#,'WEIGHT_RF']
if type[:3] == 'BGS':
    rcols.append('flux_r_dered')

if mkclusran:
    print('doing clustering randoms (possibly a 2nd time to get sys columns in)')
    tsnrcol = 'TSNR2_ELG'
    tsnrcut = 0
   
    if type[:3] == 'ELG':
        #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
        tsnrcut = 80
    if type == 'LRG':
        #dchi2 = 16  
        tsnrcut = 80  
    if type[:3] == 'BGS':
        tsnrcol = 'TSNR2_BGS'
        dchi2 = 40
        tsnrcut = 1000

    for ii in range(rm,rx):
        ct.mkclusran(dirin+type+notqso+'zdone_',dirout+type+notqso+'zdone_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits)#,ntilecut=ntile,ccut=ccut)

    

if args.nz == 'y':
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
#     if rcut is not None:
#         wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
#     if ntile > 0:
#         wzm += '_ntileg'+str(ntilecut)+'_'    
#     if ccut is not None:
#         wzm += '_'+ccut #you could change this to however you want the file names to turn out

#    regl = ['_DN','_DS','','_N','_S']
    
    if type == 'QSO':
        #zmin = 0.6
        #zmax = 4.5
        dz = 0.05
        P0 = 6000
        
    else:    
        dz = 0.02
        #zmin = 0.01
        #zmax = 1.61
    
    if type[:3] == 'LRG':
        P0 = 10000
    if type[:3] == 'ELG':
        P0 = 4000
    if type[:3] == 'BGS':
        P0 = 7000
    
    for reg in regl:
        fb = dirout+type+notqso+'zdone'+wzm+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.txt'
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0)

if args.swapz == 'y':
    import LSS.blinding_tools as blind
    for reg in regl:
        fb = dirout+type+notqso+'zdone'+reg+'_clustering.dat.fits'
        data = Table(fitsio.read(fb))
        blind.swap_z(data,fb,frac=0.01)        